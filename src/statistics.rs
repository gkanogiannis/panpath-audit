use crate::sequence::is_ambiguous;

#[derive(Clone, Copy, Default)]
/// Counts exact unit-cost edit operations and ambiguous aligned columns.
pub struct EditCounts {
    pub matched: u64,
    pub substituted: u64,
    pub source_only: u64,
    pub graph_only: u64,
    pub ambiguous_columns: u64,
}

impl EditCounts {
    pub fn edit_distance(self) -> u64 {
        self.substituted + self.source_only + self.graph_only
    }

    pub fn add(&mut self, other: Self) {
        self.matched += other.matched;
        self.substituted += other.substituted;
        self.source_only += other.source_only;
        self.graph_only += other.graph_only;
        self.ambiguous_columns += other.ambiguous_columns;
    }
}

/// Result of an exact alignment attempt under a cell budget.
pub enum Alignment {
    Complete { counts: EditCounts, cells: u64 },
    LimitExceeded,
}

#[derive(Clone, Copy)]
struct Cell {
    source: usize,
    counts: EditCounts,
}

/// Aligns source and graph sequences with a bounded wavefront computation.
pub fn align(source: &[u8], graph: &[u8], max_cells: u64) -> Alignment {
    if max_cells == 0 {
        return Alignment::LimitExceeded;
    }
    let mut cells = 1;
    let initial = extend(source, graph, 0, 0, EditCounts::default());
    if initial.source == source.len() && source.len() == graph.len() {
        return Alignment::Complete {
            counts: initial.counts,
            cells,
        };
    }
    let mut previous = vec![Some(initial)];
    let target = source.len() as i128 - graph.len() as i128;

    for score in 1..=source.len().saturating_add(graph.len()) {
        let width = score.saturating_mul(2).saturating_add(1);
        let mut current = vec![None; width];
        for k in -(score as i128)..=score as i128 {
            cells += 1;
            if cells > max_cells {
                return Alignment::LimitExceeded;
            }
            let mut best = None;
            consider(
                &mut best,
                prior(&previous, score - 1, k),
                source,
                graph,
                k,
                Operation::Substitute,
            );
            consider(
                &mut best,
                prior(&previous, score - 1, k - 1),
                source,
                graph,
                k,
                Operation::SourceOnly,
            );
            consider(
                &mut best,
                prior(&previous, score - 1, k + 1),
                source,
                graph,
                k,
                Operation::GraphOnly,
            );
            let Some(cell) = best else { continue };
            let index = usize::try_from(k + score as i128).unwrap();
            current[index] = Some(cell);
            if k == target {
                let graph_position = cell.source as i128 - k;
                if cell.source == source.len() && graph_position == graph.len() as i128 {
                    return Alignment::Complete {
                        counts: cell.counts,
                        cells,
                    };
                }
            }
        }
        previous = current;
    }
    Alignment::LimitExceeded
}

#[derive(Clone, Copy)]
enum Operation {
    Substitute,
    SourceOnly,
    GraphOnly,
}

fn prior(layer: &[Option<Cell>], score: usize, k: i128) -> Option<Cell> {
    if k < -(score as i128) || k > score as i128 {
        return None;
    }
    layer[usize::try_from(k + score as i128).ok()?]
}

fn consider(
    best: &mut Option<Cell>,
    prior: Option<Cell>,
    source: &[u8],
    graph: &[u8],
    k: i128,
    operation: Operation,
) {
    let Some(mut cell) = prior else { return };
    let graph_position = cell.source as i128
        - match operation {
            Operation::Substitute => k,
            Operation::SourceOnly => k - 1,
            Operation::GraphOnly => k + 1,
        };
    let (next_source, next_graph) = match operation {
        Operation::Substitute => (cell.source.checked_add(1), graph_position.checked_add(1)),
        Operation::SourceOnly => (cell.source.checked_add(1), Some(graph_position)),
        Operation::GraphOnly => (Some(cell.source), graph_position.checked_add(1)),
    };
    let (Some(next_source), Some(next_graph)) = (next_source, next_graph) else {
        return;
    };
    let Ok(next_graph) = usize::try_from(next_graph) else {
        return;
    };
    if next_source > source.len() || next_graph > graph.len() {
        return;
    }
    match operation {
        Operation::Substitute => {
            if next_source == 0 || next_graph == 0 {
                return;
            }
            let left = source[next_source - 1];
            let right = graph[next_graph - 1];
            if left.eq_ignore_ascii_case(&right) {
                return;
            }
            cell.counts.substituted += 1;
            cell.counts.ambiguous_columns += u64::from(is_ambiguous(left) || is_ambiguous(right));
        }
        Operation::SourceOnly => {
            if next_source == 0 {
                return;
            }
            cell.counts.source_only += 1;
            cell.counts.ambiguous_columns += u64::from(is_ambiguous(source[next_source - 1]));
        }
        Operation::GraphOnly => {
            if next_graph == 0 {
                return;
            }
            cell.counts.graph_only += 1;
            cell.counts.ambiguous_columns += u64::from(is_ambiguous(graph[next_graph - 1]));
        }
    }
    let candidate = extend(source, graph, next_source, next_graph, cell.counts);
    if best.is_none_or(|current| candidate.source > current.source) {
        *best = Some(candidate);
    }
}

fn extend(
    source: &[u8],
    graph: &[u8],
    mut source_position: usize,
    mut graph_position: usize,
    mut counts: EditCounts,
) -> Cell {
    while source_position < source.len()
        && graph_position < graph.len()
        && source[source_position].eq_ignore_ascii_case(&graph[graph_position])
    {
        counts.matched += 1;
        counts.ambiguous_columns +=
            u64::from(is_ambiguous(source[source_position]) || is_ambiguous(graph[graph_position]));
        source_position += 1;
        graph_position += 1;
    }
    Cell {
        source: source_position,
        counts,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn complete(source: &str, graph: &str) -> EditCounts {
        match align(source.as_bytes(), graph.as_bytes(), 10_000) {
            Alignment::Complete { counts, .. } => counts,
            Alignment::LimitExceeded => panic!("alignment limit"),
        }
    }

    #[test]
    fn counts_edit_operations() {
        let substitution = complete("ACGT", "ACTT");
        assert_eq!((substitution.matched, substitution.substituted), (3, 1));

        let source_only = complete("ACGT", "AGT");
        assert_eq!((source_only.matched, source_only.source_only), (3, 1));

        let graph_only = complete("AGT", "ACGT");
        assert_eq!((graph_only.matched, graph_only.graph_only), (3, 1));
    }

    #[test]
    fn counts_literal_ambiguous_columns() {
        let counts = complete("AN", "AN");
        assert_eq!(counts.matched, 2);
        assert_eq!(counts.ambiguous_columns, 1);
    }

    #[test]
    fn enforces_cell_limit() {
        assert!(matches!(
            align(b"AAAA", b"TTTT", 1),
            Alignment::LimitExceeded
        ));
    }

    #[test]
    fn agrees_with_dynamic_programming_on_small_inputs() {
        let values = ["", "A", "T", "AA", "AT", "TA", "TT", "AAT", "TTA"];
        for source in values {
            for graph in values {
                let counts = complete(source, graph);
                assert_eq!(
                    counts.edit_distance(),
                    distance(source.as_bytes(), graph.as_bytes()),
                    "{source:?} {graph:?}"
                );
                assert_eq!(
                    counts.matched + counts.substituted + counts.source_only,
                    source.len() as u64
                );
                assert_eq!(
                    counts.matched + counts.substituted + counts.graph_only,
                    graph.len() as u64
                );
            }
        }
    }

    fn distance(left: &[u8], right: &[u8]) -> u64 {
        let mut previous: Vec<u64> = (0..=right.len() as u64).collect();
        for (i, left_base) in left.iter().enumerate() {
            let mut current = vec![i as u64 + 1];
            for (j, right_base) in right.iter().enumerate() {
                current.push(
                    (previous[j + 1] + 1)
                        .min(current[j] + 1)
                        .min(previous[j] + u64::from(!left_base.eq_ignore_ascii_case(right_base))),
                );
            }
            previous = current;
        }
        previous[right.len()]
    }
}
