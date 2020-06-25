const width = 1500;
const height = 1500;
const MAX_DISTANCE = 3;
const MIN_DISTANCE = 0.1;
const SCALE = 10;
const INNER_RADIUS = 300;
const RADIUS_TICK = 50;
const svg = d3
  .select(".figure")
  .append("svg")
  .attr("width", width)
  .attr("height", height);
const graph_layer = svg
  .append("g")
  .attr("transform", `translate(${height / 2},${width / 2})`);

const info = d3.select(".info");

const drag = (simulation) => {
  function dragstarted(d) {
    if (!d3.event.active) simulation.alphaTarget(0.3).restart();
    d.fx = d.x;
    d.fy = d.y;
  }
  function dragged(d) {
    d.fx = d3.event.x;
    d.fy = d3.event.y;
  }
  function dragended(d) {
    if (!d3.event.active) simulation.alphaTarget(0);
    d.fx = null;
    d.fy = null;
  }
  return d3
    .drag()
    .on("start", dragstarted)
    .on("drag", dragged)
    .on("end", dragended);
};

const plotData = (dataset, cluster_num) =>
  d3
    .json(dataset)
    .then((dataset) => {
      const [units, eds] = dataset;
      const max_unit = Math.max(...units.map((d) => d.id));
      const edges = eds
        .map((d) => {
          d[0].count = d[1];
          d[0].source = d[0].from + d[0].from_cluster * max_unit;
          d[0].target = d[0].to + d[0].to_cluster * max_unit;
          return d[0];
        })
        .filter((d) => d.count > 3);
      const radius_scale = d3
        .scaleLinear()
        .domain([0, max_unit])
        .range([0, 2 * Math.PI]);
      const useful_node = new Set(edges.flatMap((d) => [d.source, d.target]));
      const units_with_cluster = Array.from(
        { length: cluster_num },
        (_, i) => i
      )
        .flatMap((i) =>
          units.map((d) => {
            const theta = radius_scale(d.position);
            const r = INNER_RADIUS + d.cluster * RADIUS_TICK;
            return {
              id: d.id,
              position: d.position,
              cluster: i,
              unit_id: d.id + max_unit * i,
              x: Math.cos(theta) * r,
              y: Math.sin(theta) * r,
            };
          })
        )
        .filter((d) => useful_node.has(d.unit_id));
      const node_layer = graph_layer
        .append("g")
        .attr("class", "nodes")
        .selectAll(".node")
        .data(units_with_cluster)
        .enter()
        .append("circle")
        .attr("class", "node")
        .attr("fill", "black")
        .attr("r", 4);
      //     .attr("cx", (d) => {
      //   const theta = radius_scale(d.position);
      //   const r = INNER_RADIUS + d.cluster * RADIUS_TICK;
      //   return Math.cos(theta) * r;
      // })
      // .attr("cy", (d) => {
      //   const theta = radius_scale(d.position);
      //   const r = INNER_RADIUS + d.cluster * RADIUS_TICK;
      //   return Math.sin(theta) * r;
      // });
      const max_width = Math.max(...edges.map((d) => d.count));
      const stroke_scale = d3
        .scaleLinear()
        .domain([0, max_width])
        .range([0, 10]);
      const id_to_position = new Map(units.map((d) => [d.id, d.position]));
      const link_layer = graph_layer
        .append("g")
        .attr("class", "links")
        .selectAll(".link")
        .data(edges)
        .enter()
        .append("line")
        .attr("class", "link")
        .attr("stroke", (d) => d3.schemeSet1[d.answer])
        .attr("stroke-width", (d) => stroke_scale(d.count));
      // .attr("x1", (d) => {
      //   const theta = radius_scale(id_to_position.get(d.from));
      //   const r = INNER_RADIUS + d.from_cluster * RADIUS_TICK;
      //   return Math.cos(theta) * r;
      // })
      // .attr("y1", (d) => {
      //   const theta = radius_scale(id_to_position.get(d.from));
      //   const r = INNER_RADIUS + d.from_cluster * RADIUS_TICK;
      //   return Math.sin(theta) * r;
      // })
      // .attr("x2", (d) => {
      //   const theta = radius_scale(id_to_position.get(d.to));
      //   const r = INNER_RADIUS + d.to_cluster * RADIUS_TICK;
      //   return Math.cos(theta) * r;
      // })
      // .attr("y2", (d) => {
      //   const theta = radius_scale(id_to_position.get(d.to));
      //   const r = INNER_RADIUS + d.to_cluster * RADIUS_TICK;
      //   return Math.sin(theta) * r;
      // });
      const distance_scale = d3
        .scaleLinear()
        .domain([0, max_width])
        .range([MAX_DISTANCE, MIN_DISTANCE]);
      const simulation = d3
        .forceSimulation()
        .alphaDecay(0.001)
        .nodes(units_with_cluster)
        .force(
          "charge",
          d3.forceManyBody().strength(() => -0.06)
        )
        .force(
          "link",
          d3
            .forceLink()
            .links(edges)
            .id((d) => d.unit_id)
            .distance((d) => distance_scale(d.count))
          //.strength((d) => power_scale(d.count))
        )
        .force("center", d3.forceCenter())
        .on("tick", () => {
          link_layer
            .attr("x1", (l) => l.source.x)
            .attr("y1", (l) => l.source.y)
            .attr("x2", (l) => l.target.x)
            .attr("y2", (l) => l.target.y);
          node_layer.attr("cx", (d) => d.x).attr("cy", (d) => d.y);
        });
      node_layer.call(drag(simulation));
    })
    .then(
      (ok) => ok,
      (why) => console.log(why)
    );
