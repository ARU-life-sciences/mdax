const data = IRX_DATA.map((d) => ({
  ...d,
  tir_ident: d.tir_ident ?? NaN
}));

function formatBp(x) {
  if (x == null || !isFinite(x)) return "";
  const ax = Math.abs(x);
  if (ax >= 1e9) return `${d3.format(".2f")(x / 1e9)}Gb`;
  if (ax >= 1e6) return `${d3.format(".2f")(x / 1e6)}Mb`;
  if (ax >= 1e3) return `${d3.format(".1f")(x / 1e3)}Kb`;
  return `${x} bp`;
}

function containerWidth(id) {
  return document.getElementById(id).clientWidth;
}

const formatInt = d3.format(",");

const state = {
  minIdent: 0.5,
  minKept: 50,
  minMatches: 50,
  minArm: 500,
  minSpacer: 0,
  classes: new Set(["immediate", "spaced", "wide", "overlap"])
};

function filteredData() {
  return data.filter((d) =>
    (isFinite(d.tir_ident) ? d.tir_ident : -Infinity) >= state.minIdent &&
    d.kept_pts >= state.minKept &&
    d.matches >= state.minMatches &&
    d.arm_len >= state.minArm &&
    d.spacer >= state.minSpacer &&
    state.classes.has(d.ir_class)
  );
}

function renderControls() {
  const root = document.getElementById("controls");
  root.innerHTML = "";

  const controls = [
    ["minIdent", "min tir_ident", 0, 1, 0.01],
    ["minKept", "min kept_pts", 0, 2000, 1],
    ["minMatches", "min matches", 0, 5000, 1],
    ["minArm", "min arm_len (bp)", 0, 200000, 10],
    ["minSpacer", "min spacer (bp)", 0, 200000, 10]
  ];

  for (const [key, label, min, max, step] of controls) {
    const wrap = document.createElement("div");
    const lab = document.createElement("label");
    const val = document.createElement("div");
    const input = document.createElement("input");

    lab.textContent = label;
    val.textContent = state[key];
    input.type = "range";
    input.min = min;
    input.max = max;
    input.step = step;
    input.value = state[key];
    input.addEventListener("input", () => {
      state[key] = +input.value;
      val.textContent = state[key];
      renderAll();
    });

    wrap.appendChild(lab);
    wrap.appendChild(input);
    wrap.appendChild(val);
    root.appendChild(wrap);
  }

  const classWrap = document.createElement("div");
  classWrap.innerHTML = `<label>ir_class</label>`;
  for (const cls of ["overlap", "immediate", "spaced", "wide", "unknown"]) {
    const row = document.createElement("div");
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.checked = state.classes.has(cls);
    cb.addEventListener("change", () => {
      if (cb.checked) state.classes.add(cls);
      else state.classes.delete(cls);
      renderAll();
    });
    row.appendChild(cb);
    row.append(` ${cls}`);
    classWrap.appendChild(row);
  }
  root.appendChild(classWrap);
}

function renderPlot(id, plot) {
  const node = document.getElementById(id);
  node.innerHTML = "";
  node.appendChild(plot);
}

function renderAll() {
  const filtered = filteredData();

  // Architecture plot
  renderPlot("plot-arch", Plot.plot({
    width: containerWidth("plot-arch"),
    height: 520,
    marginLeft: 65,
    marginBottom: 50,
    x: { label: "Spacer (bp)", grid: true, type: "log", labelOffset: 40 },
    y: { label: "Arm length (bp)", grid: true, type: "log", labelOffset: 0 },
    color: { label: "tir_ident", legend: true },
    marks: [
      Plot.dot(filtered.filter((d) => d.spacer > 0 && d.arm_len > 0), {
        x: "spacer",
        y: "arm_len",
        r: (d) => Math.sqrt(Math.max(1, d.kept_pts)) / 4,
        fill: "tir_ident",
        fillOpacity: 0.85
      }),
      Plot.tip(
        filtered.filter((d) => d.spacer > 0 && d.arm_len > 0),
        Plot.pointer({
          x: "spacer",
          y: "arm_len",
          title: (d) => [
            `Contig: ${d.contig}`,
            `Interval: ${formatBp(d.start)} – ${formatBp(d.end)}`,
            `Arm length: ${formatBp(d.arm_len)}`,
            `Spacer: ${formatBp(d.spacer)}`,
            `Left arm: ${formatBp(d.left_arm_len)}`,
            `Right arm: ${formatBp(d.right_arm_len)}`,
            `TIR identity: ${d.tir_ident.toFixed(3)}`,
            `kept_pts: ${formatInt(d.kept_pts)}`,
            `matches: ${formatInt(d.matches)}`,
            `class: ${d.ir_class}`
          ].join("\n")
        })
      )
    ]
  }));

  // Arm symmetry plot
  const symData = filtered.filter(
    (d) => d.left_arm_len > 0 && d.right_arm_len > 0
  );

  const armMax = d3.max(symData, d => Math.max(d.left_arm_len, d.right_arm_len));
  const armMin = Math.max(
    1,
    d3.min(symData, d => Math.min(d.left_arm_len, d.right_arm_len))
  );

  renderPlot("plot-sym", Plot.plot({
    width: containerWidth("plot-sym"),
    height: 440,
    marginBottom: 50,
    x: {
      label: "Left arm length (bp)",
      grid: true,
      type: "log",
      labelOffset: 40,
      domain: [armMin, armMax],
    },
    y: {
      label: "Right arm length (bp)",
      grid: true,
      type: "log",
      labelOffset: 0,
      domain: [armMin, armMax],
    },
    color: {
      label: "Inverted repeat % identity",
      legend: true
    },
    marks: [
      Plot.line(
        [
          [armMin, armMin],
          [armMax, armMax]
        ],
        {
          x: (d) => d[0],
          y: (d) => d[1],
          strokeOpacity: 0.3,
          strokeDasharray: "4,4"
        }
      ),
      Plot.dot(
        filtered.filter((d) => d.left_arm_len > 0 && d.right_arm_len > 0),
        {
          x: "left_arm_len",
          y: "right_arm_len",
          fill: "tir_ident",
          fillOpacity: 0.8,
          r: (d) => Math.max(2, Math.sqrt(d.kept_pts) / 5)
        }
      ),
      Plot.tip(
        filtered.filter((d) => d.spacer > 0 && d.arm_len > 0),
        Plot.pointer({
          x: "left_arm_len",
          y: "right_arm_len",

          title: (d) =>
            [
              `Contig: ${d.contig}`,
              `Interval: ${formatBp(d.start)} – ${formatBp(d.end)}`,
              `Arm length: ${formatBp(d.arm_len)}`,
              `Spacer: ${formatBp(d.spacer)}`,
              `Left arm: ${formatBp(d.left_arm_len)}`,
              `Right arm: ${formatBp(d.right_arm_len)}`,
              `TIR identity: ${d.tir_ident.toFixed(3)}`,
              `kept_pts: ${formatInt(d.kept_pts)}`,
              `matches: ${formatInt(d.matches)}`,
              `class: ${d.ir_class}`
            ].join("\n")
        })
      )
    ]
  }))

  // Identity distribution
  renderPlot("plot-ident", Plot.plot({
    width: containerWidth("plot-ident"),
    height: 260,
    marginBottom: 50,
    x: { label: "Inverted repeat % identity", grid: true, labelOffset: 40 },
    y: { label: "Count", labelOffset: 0 },
    color: { legend: true },
    marks: [
      Plot.rectY(
        filtered,
        Plot.binX(
          { y: "count" },
          { x: "tir_ident", fill: "ir_class", thresholds: 30 }
        )
      ),
      Plot.ruleY([0])
    ]
  }))

  // 
  renderPlot("plot-interval", Plot.plot({
    width: containerWidth("plot-interval"),
    height: 260,
    marginBottom: 50,
    x: {
      label: "Interval length (bp)",
      grid: true,
      type: "log",
      labelOffset: 40
    },
    y: { label: "Count", labelOffset: 0 },
    color: { legend: true },
    marks: [
      Plot.rectY(
        filtered.filter((d) => d.interval_len > 0),
        Plot.binX(
          { y: "count" },
          { x: "interval_len", fill: "ir_class", thresholds: 30 }
        )
      ),
      Plot.ruleY([0])
    ]
  }))

  // Plot kept
  renderPlot("plot-kept", Plot.plot({
    width: containerWidth("plot-kept"),
    height: 420,
    marginBottom: 50,
    x: { label: "Inverted repeat % identity", grid: true, labelOffset: 40 },
    y: { label: "Kept points", grid: true, type: "log", labelOffset: 0 },
    color: { legend: true },
    marks: [
      Plot.dot(filtered, {
        x: "tir_ident",
        y: "kept_pts",
        fill: "ir_class",
        r: (d) => Math.max(2, Math.sqrt(d.arm_len) / 20),
        fillOpacity: 0.8
      }),
      Plot.tip(
        filtered,
        Plot.pointer({
          x: "tir_ident",
          y: "kept_pts",

          title: (d) =>
            [
              `Contig: ${d.contig}`,
              `Break: ${formatBp(d.break_pos)}`,
              `TIR identity: ${d.tir_ident.toFixed(3)}`,
              `Arm length: ${formatBp(d.arm_len)}`,
              `Spacer: ${formatBp(d.spacer)}`,
              `kept_pts: ${formatInt(d.kept_pts)}`,
              `matches: ${formatInt(d.matches)}`,
              `class: ${d.ir_class}`
            ].join("\n")
        })
      )
    ]
  }))

  // Plot matches
  renderPlot("plot-matches", Plot.plot({
    width: containerWidth("plot-matches"),
    height: 420,
    marginBottom: 50,
    x: { label: "Inverted repeat % identity", grid: true, labelOffset: 40 },
    y: { label: "Matches", grid: true, type: "log", labelOffset: 0 },
    color: { legend: true },
    marks: [
      Plot.dot(filtered, {
        x: "tir_ident",
        y: "matches",
        fill: "ir_class",
        fillOpacity: 0.8
      }),

      Plot.tip(
        filtered,
        Plot.pointer({
          x: "tir_ident",
          y: "matches",

          title: (d) =>
            [
              `Contig: ${d.contig}`,
              `Break: ${formatBp(d.break_pos)}`,
              `TIR identity: ${d.tir_ident.toFixed(3)}`,
              `Arm length: ${formatBp(d.arm_len)}`,
              `Spacer: ${formatBp(d.spacer)}`,
              `kept_pts: ${formatInt(d.kept_pts)}`,
              `matches: ${formatInt(d.matches)}`,
              `class: ${d.ir_class}`
            ].join("\n")
        })
      )
    ]
  }))

  // plot the IR's on the contigs
  const topContigs = Array.from(
    d3.rollup(
      filtered,
      (v) => ({
        n: v.length,
        contig_len: v[0].contig_len
      }),
      (d) => d.contig
    ),
    ([contig, stats]) => ({ contig, ...stats })
  ).sort((a, b) => d3.descending(a.contig_len, b.contig_len));

  const contigSet = new Set(topContigs.map((d) => d.contig));

  renderPlot("plot-contigs", Plot.plot({
    width: containerWidth("plot-contigs"),
    height: 52 * Math.max(6, topContigs.length),
    style: "overflow: visible;",
    marginLeft: 120,
    marginBottom: 50,
    x: { label: "Contig coordinate (bp)", grid: true, tickFormat: formatBp, labelOffset: 40, domain: [0, d3.max(topContigs, d => d.contig_len)] },
    y: {
      label: null,
      domain: topContigs.map((d) => d.contig),
      tickSize: 0,
    },
    color: {
      label: "Inverted repeat % identity",
      legend: true
    },
    marks: [
      Plot.ruleY(
        topContigs.map((d) => d.contig),
        { strokeOpacity: 0.12 }
      ),
      Plot.ruleX([0], { strokeOpacity: 0.3 }),

      // contig end marker
      Plot.tickX(topContigs, {
        x: "contig_len",
        y: "contig",
        stroke: "#555",
        strokeWidth: 2
      }),

      // IR intervals
      Plot.ruleX(
        filtered.filter((d) => contigSet.has(d.contig)),
        {
          x1: "start",
          x2: "end",
          y: "contig",
          stroke: "tir_ident",
          strokeWidth: 2.5
        }
      ),

      // breakpoints
      Plot.tickX(filtered, {
        x: "break_pos",
        y: "contig",
        stroke: "tir_ident",
        strokeOpacity: 0.45
      }),

      Plot.tip(
        topContigs,
        Plot.pointer({
          x: "contig_len",
          y: "contig",
          title: (d) =>
            [
              `Contig: ${d.contig}`,
              `Contig length: ${formatBp(d.contig_len)}`,
              `IR count: ${formatInt(d.n)}`
            ].join("\n")
        })
      ),

      Plot.tip(
        filtered,
        Plot.pointer({
          x: "break_pos",
          y: "contig",
          title: (d) =>
            [
              `Contig: ${d.contig}`,
              `Contig length: ${formatBp(d.contig_len)}`,
              `Interval: ${formatBp(d.start)} – ${formatBp(d.end)}`,
              `Length: ${formatBp(d.interval_len)}`,
              `Breakpoint: ${formatBp(d.break_pos)}`,
              `Arm length: ${formatBp(d.arm_len)}`,
              `Spacer: ${formatBp(d.spacer)}`,
              `TIR identity: ${d.tir_ident.toFixed(3)}`,
              `kept_pts: ${formatInt(d.kept_pts)}`,
              `matches: ${formatInt(d.matches)}`,
              `class: ${d.ir_class}`
            ].join("\n")
        })
      )
    ]
  })
  )

}

renderControls();
renderAll();

window.addEventListener("resize", () => {
  renderAll();
});
