const width = 2000;
const height = 1000;
const svg = d3
      .select("#plot")
      .append("svg")
      .attr("width",width)
      .attr("height", height)
      .attr("viewBox", [-width/2,-height/2,width,height]);

const graph_layer = svg
      .append("g")
      .attr("class", "hap-network");

const link_svg = svg
      .append("g")
      .attr("stroke","black")
      .attr("stroke-opacity", 0.5);

const node_svg = svg
      .append("g")
      .attr("fill", "gray")
      .attr("stroke","#999")
      .attr("stroke-width",1.5);
const CIRCLE_SIZE = 6;
const EDGE_WEIGHT = 4;


const drag = simulation => {
    const dragStarted = (event, d) => {
        if (!event.active) simulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
    };
    const dragged = (event, d) => {
        d.fx = event.x;
        d.fy = event.y;
    };
    const dragEnded = (event, d) => {
        if (!event.active) simulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
    };
    return d3.drag()
        .on("start", dragStarted)
        .on("drag", dragged)
        .on("end", dragEnded);
};

const plotData = (path_to_json)=>
      Promise.all([path_to_json].map((file) => d3.json(file)))
      .then(([graph]) => {
          const getFraction = (node)=> node.maternal / (node.maternal + node.parental);
          const nodes = graph.nodes;
          const inNodes = new Set(nodes.map(n => n.id));
          const edges = graph.edges.filter(edge=> edge.source <= edge.target)
                .filter(edge =>  inNodes.has(edge.source) && inNodes.has(edge.target));
          const maxCoverage = Math.max(...nodes.map(n=>n.maternal + n.parental));
          const maxWeight = Math.max(...edges.map(e=>Math.abs(e.weight)));
          const target_w = 200;
          const forceLink = d3.forceLink(edges)
                .id((d)=>d.id)
                .strength(d=> d.weight/maxWeight);
          const simulation = d3.forceSimulation(nodes)
                .force("link", forceLink)
                .force("charge",d3.forceManyBody().strength(-10))
                .force("x",d3.forceX().strength(0.05))
                .force("y", d3.forceY());
                //.force("y", d3.forceY().y(d => (getFraction(d)-0.5) * target_w));
          console.log(`nodes:${nodes.length}`);
          console.log(`edges:${edges.length}`);
          const link = link_svg.selectAll("line")
                .data(edges)
                .join("line")
                .attr("stroke-width", d=> Math.abs(d.weight)/maxWeight * EDGE_WEIGHT)
                .attr("stroke", d=> d.weight > 0 ? "gray" : "blue");
                // .attr("visibility", d=> d.weight > 0 ? "visible" : "hidden");
          const node = node_svg.selectAll("circle")
                .data(nodes)
                .join("circle")
                .attr("r",d=>(d.maternal + d.parental)/maxCoverage * CIRCLE_SIZE)
          //.attr("fill", d => d.maternal < d.parental ? "red" : "blue");
                .attr("fill", d => d3.interpolateViridis(getFraction(d)))
                .call(drag(simulation));
          
          simulation.on("tick",()=>{
              link
                  .attr("x1", d => d.source.x)
                  .attr("y1",d => d.source.y)
                  .attr("x2",d => d.target.x)
                  .attr("y2", d => d.target.y);
              node
                  .attr("cx", d => d.x)
                  .attr("cy", d => d.y);
          });
          return svg.nodes();
      })
      .then(
          (ok) => ok,
          (why)=>console.log(why)
      );


