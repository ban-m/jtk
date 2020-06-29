fn main() {
    let nodes = 12;
    let edges = vec![
        vec![1, 5, 6, 9],
        vec![0, 2, 4],
        vec![1, 4, 3],
        vec![2, 7, 8, 10],
        vec![1, 2, 6, 7],
        vec![0, 9],
        vec![0, 4, 9],
        vec![3, 4, 8, 10],
        vec![3, 7, 10, 11],
        vec![0, 5, 6],
        vec![3, 7, 8, 11],
        vec![8, 10],
    ];
    let edges: Vec<Vec<_>> = edges
        .into_iter()
        .map(|eds| eds.into_iter().map(|x| (x, 1)).collect())
        .collect();
    let graph = mcl::Graph::new(nodes, &edges);
    eprintln!("Result:{:?}", graph.clustering(2, 2));
}
