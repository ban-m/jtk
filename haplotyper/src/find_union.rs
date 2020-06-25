#![allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct FindUnion {
    /// The vector of parents. If parents[i] = j, the j-th node represnets
    /// the component which i-th node redides.
    parents: Vec<usize>,
    /// The size of each component. If sizes[i] = k, the size of component which i-th nodes resides
    /// is k.
    /// It is valid if and only if i-th nodes is the represnetor of the components.
    sizes: Vec<Option<usize>>,
    /// Length of the dataset.
    length: usize,
}

impl FindUnion {
    /// Create a new instance of FindUnion tree.
    pub fn new(size: usize) -> Self {
        let parents = (0..size).collect();
        let sizes = (0..size).map(|_| Some(1)).collect();
        let length = size;
        FindUnion {
            parents,
            sizes,
            length,
        }
    }
    /// Clear the contents of the instance.
    pub fn clear(&mut self) {
        self.parents.clear();
        self.sizes.clear();
        self.length = 0;
    }
    /// Check if the contents is empty.
    pub fn is_empty(&self) -> bool {
        (self.length == 0) && self.parents.is_empty() && self.sizes.is_empty()
    }
    /// Recreate the instance. Note that it would erase all the exiting contents of this instance.
    pub fn refresh(&mut self, size: usize) {
        self.clear();
        assert!(self.is_empty());
        self.length = size;
        for i in 0..size {
            self.parents.push(i);
            self.sizes.push(Some(1))
        }
    }
    /// Find the representative of nodes [index].
    /// If the argument goes out of range, return None,
    pub fn find(&mut self, index: usize) -> Option<usize> {
        if index > self.length {
            None
        } else {
            let parent = {
                let mut temp = index;
                while temp != self.parents[temp] {
                    temp = self.parents[temp];
                }
                temp
            };
            let mut index = index;
            while index != parent {
                let next = self.parents[index];
                self.parents[index] = parent;
                index = next;
            }
            Some(parent)
        }
    }
    /// Unite the comporents node1 residing and the components
    /// nodes2 residing.
    /// If either of index exceeds the range, return None,
    /// else and if the marging process success, return Some(true).
    /// else, panic through unreachable!().
    pub fn unite(&mut self, node1: usize, node2: usize) -> Option<()> {
        if node1 >= self.length || node2 >= self.length {
            None
        } else if node1 == node2 {
            // Need not to merge.
            Some(())
        } else {
            let parent1 = self.find(node1)?;
            let parent2 = self.find(node2)?;
            if parent1 != parent2 {
                if self.sizes[parent1] > self.sizes[parent2] {
                    self.parents[parent2] = parent1;
                    self.sizes[parent1] = Some(self.sizes[parent1]? + self.sizes[parent2]?);
                    self.sizes[parent2] = None;
                } else {
                    self.parents[parent1] = parent2;
                    self.sizes[parent2] = Some(self.sizes[parent1]? + self.sizes[parent2]?);
                    self.sizes[parent1] = None;
                }
            }
            Some(())
        }
    }
    #[allow(dead_code)]
    /// Determine if the node1 and node2 resides in the same component.
    /// return None if either of index exceeds the range,
    pub fn same(&mut self, node1: usize, node2: usize) -> Option<bool> {
        if node1 >= self.length || node2 >= self.length {
            None
        } else {
            Some(self.find(node1) == self.find(node2))
        }
    }
    /// Get the size of the component nodes1 resides.
    /// Return None when out of range.
    pub fn size(&mut self, node1: usize) -> Option<usize> {
        if node1 >= self.length {
            None
        } else {
            let parent = self.find(node1).unwrap();
            self.sizes[parent]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn findunion_init() {
        FindUnion::new(0);
    }

    #[test]
    fn findunion_init2() {
        FindUnion::new(10);
    }

    #[test]
    fn unite() {
        let mut fu = FindUnion::new(10);
        assert_eq!(Some(()), fu.unite(1, 3));
        assert_eq!(None, fu.unite(100, 0));
        assert_eq!(None, fu.unite(10, 3));
        assert_eq!(None, fu.unite(3, 10));
    }

    #[test]
    fn find() {
        let mut fu = FindUnion::new(10);
        fu.unite(0, 1);
        assert_eq!(fu.find(0), fu.find(1));
        fu.unite(2, 3);
        assert_eq!(fu.find(2), fu.find(3));
    }

    #[test]
    fn find2() {
        let mut fu = FindUnion::new(10);
        fu.unite(0, 1);
        fu.unite(1, 2);
        assert_eq!(fu.find(0), fu.find(2));
    }

    #[test]
    fn find3() {
        let mut fu = FindUnion::new(10);
        fu.unite(0, 1);
        fu.unite(0, 2);
        assert_eq!(fu.find(0), fu.find(2));
    }

    #[test]
    fn find4() {
        let mut fu = FindUnion::new(10);
        fu.unite(0, 1);
        fu.unite(4, 2);
        fu.unite(1, 4);
        assert_eq!(fu.find(0), fu.find(2));
    }
    #[test]
    fn same() {
        let mut fu = FindUnion::new(10);
        fu.unite(0, 1);
        assert_eq!(fu.same(0, 1), Some(true));
        fu.unite(2, 3);
        assert_eq!(fu.same(2, 3), Some(true));
        assert_eq!(fu.same(0, 8), Some(false));
        assert_eq!(fu.same(213, 232), None);
    }

    #[test]
    fn same2() {
        let mut fu = FindUnion::new(10);
        fu.unite(0, 1);
        fu.unite(1, 2);
        assert_eq!(fu.same(0, 2), Some(true));
    }

    #[test]
    fn same3() {
        let mut fu = FindUnion::new(10);
        fu.unite(0, 1);
        fu.unite(0, 2);
        assert_eq!(fu.same(0, 2), Some(true));
    }

    #[test]
    fn same4() {
        let mut fu = FindUnion::new(10);
        fu.unite(0, 1);
        fu.unite(4, 2);
        fu.unite(1, 4);
        assert_eq!(fu.same(0, 2), Some(true));
    }

    #[test]
    fn size() {
        let mut fu = FindUnion::new(10);
        for i in 0..10 {
            assert_eq!(fu.size(i), Some(1));
        }
        assert_eq!(fu.size(100), None);
    }

    #[test]
    fn size2() {
        let mut fu = FindUnion::new(10);
        fu.unite(0, 1);
        assert_eq!(fu.size(1), Some(2));
        assert_eq!(fu.size(0), Some(2));
        assert_eq!(fu.size(3), Some(1));
    }

    #[test]
    fn size3() {
        let mut fu = FindUnion::new(10);
        fu.unite(0, 1);
        fu.unite(0, 2);
        fu.unite(4, 3);
        fu.unite(3, 2);
        for i in 0..5 {
            assert_eq!(fu.size(i), Some(5));
        }
    }

    #[test]
    fn size4() {
        let mut fu = FindUnion::new(10);
        for i in 0..9 {
            fu.unite(i, i + 1);
        }
        for i in 0..10 {
            assert_eq!(fu.size(i), Some(10));
        }
    }
    #[test]
    fn all_connected() {
        let mut fu = FindUnion::new(10);
        for i in 0..10 {
            for j in (0..10).rev() {
                fu.unite(i, j);
                fu.unite((i * 3) % 10, (j * 7) % 10);
            }
        }
        let p = fu.find(0).unwrap();
        assert!((0..10).all(|e| fu.find(e).unwrap() == p));
        let (_, size) = (0..10)
            .filter_map(|e| match fu.find(e) {
                Some(p) if p == e => Some((e, fu.size(e)?)),
                _ => None,
            })
            .max_by_key(|&(_, size)| size)
            .unwrap();
        assert_eq!(size, 10);
    }
}
