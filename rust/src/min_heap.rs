/// A min-heap, this differs from rusts BinaryHeap
/// in that we can remove any item in the heap out-of-order and duplicates.
///
/// Characteristics:
///
/// - Uses {value: user_data}, where the 'value' must support PartialOrd
///   for ordering in the heap.
/// - Supported duplicate entries,
///   (Note that the order, while not *undefined* is determined by the binary tree structure).
///
/// Overview:
///
/// Module:
/// - MinHeap::new() -> MinHeap
/// - MinHeap::with_capacity(capacity) -> MinHeap
///
/// Methods:
/// - heap.insert(sort_value, user_data) -> handle
/// - heap.insert_or_update(handle, sort_value, user_data)
/// - heap.node_value_update(handle, sort_value) -> handle
/// - heap.node_value_update_with_data(handle, sort_value, user_data) -> handle

/// - heap.remove(handle)
/// - heap.pop_min() -> Option(user_data)
///

/// Invalid index.
const INVALID: usize = ::std::usize::MAX;

/// Use only for: `self.nodes[NodeHandle]`
/// While this is just an index internally `NodeHandle` is opaque
/// to prevent external users mixing with other types.

// even though we don't want users of this struct to meddle with its internals
// its useful to be able to compare them.
#[derive(PartialEq, Debug, Copy, Clone)]
pub struct NodeHandle(usize);

impl NodeHandle {
    pub const INVALID: NodeHandle = NodeHandle(INVALID);
}

pub trait HeapValue: PartialOrd + PartialEq + Copy {}
impl<TOrd> HeapValue for TOrd where TOrd: PartialOrd + PartialEq + Copy {}

pub trait HeapData: Copy {}
impl<TData> HeapData for TData where TData: Copy {}

pub struct Node<TOrd: HeapValue, TData: HeapData> {
    /// Value to order by.
    value: TOrd,

    /// Data supplied by the API user.
    user_data: TData,

    /// index into `MinHeap.tree_index`
    ///
    /// When free'd doubles as a single-linked list into nodes,
    /// so we can re-use them.
    index: usize,
}

pub struct MinHeap<TOrd: HeapValue, TData: HeapData> {
    /// Index into `node` array.
    tree_index: Vec<usize>,

    /// Node storage, unused nodes are referenced by `free`.
    node: Vec<Node<TOrd, TData>>,

    /// Chain free nodes, where each nodes index points to the next free node
    /// terminating once an `INVALID` value is reached.
    free: usize,
}

fn bin_parent(i: usize) -> usize {
    (i - 1) >> 1
}
fn bin_left(i: usize) -> usize {
    (i << 1) + 1
}
fn bin_right(i: usize) -> usize {
    (i << 1) + 2
}

macro_rules! unlikely { ($body:expr) => { $body } }

impl<TOrd: HeapValue, TData: HeapData> MinHeap<TOrd, TData> {

    // -------------------------------------------------------------------
    // Private API
    //
    /// Compare by value, using index as tie-breaker for deterministic ordering.
    /// Otherwise equal values would have arbitrary ordering based on the heap's internal state.
    /// While technically deterministic, this is less stable as minor changes to the heap contents
    /// can yield different results elsewhere in the heap.
    fn node_compare(
        a: &Node<TOrd, TData>,
        b: &Node<TOrd, TData>,
    ) -> bool {
        (a.value < b.value) || ((a.value == b.value) && (a.index < b.index))
    }

    // Debug only, does full search on data!
    // ensures we don't allow incorrect insertion/removal.
    fn contains_node_handle(
        &self, nhandle: &NodeHandle,
    ) -> bool {
        self.tree_index.iter().any(|&i| i == nhandle.0)
    }

    /// `self.tree(i)`, short for `self.node[self.tree_index[i]]`
    #[inline(always)]
    fn tree(
        &self, i: usize,
    ) -> &Node<TOrd, TData> {
        debug_assert!(i < self.tree_index.len());
        unsafe { self.node.get_unchecked(*self.tree_index.get_unchecked(i)) }
    }
    #[allow(dead_code)]
    #[inline(always)]
    fn tree_mut(
        &mut self, i: usize,
    ) -> &mut Node<TOrd, TData> {
        debug_assert!(i < self.tree_index.len());
        unsafe { self.node.get_unchecked_mut(*self.tree_index.get_unchecked(i)) }
    }

    fn heap_swap(&mut self, i: usize, j: usize) {
        self.tree_index.swap(i, j);

        unsafe {
            let i_node = *self.tree_index.get_unchecked(i);
            let j_node = *self.tree_index.get_unchecked(j);
            let t = self.node.get_unchecked(i_node).index;
            self.node.get_unchecked_mut(i_node).index = self.node.get_unchecked(j_node).index;
            self.node.get_unchecked_mut(j_node).index = t;
        }
    }

    fn heap_compare(
        &self, i: usize, j: usize,
    ) -> bool {
        MinHeap::node_compare(self.tree(i), self.tree(j))
    }

    fn heap_down(&mut self, mut i: usize) {
        // size won't change in the loop
        let size = self.tree_index.len();

        loop {
            let l = bin_left(i);
            let r = bin_right(i);
            let mut smallest = i;

            if (l < size) && self.heap_compare(l, smallest) {
                smallest = l
            }
            if (r < size) && self.heap_compare(r, smallest) {
                smallest = r;
            }

            if smallest == i {
                break;
            }

            self.heap_swap(i, smallest);

            i = smallest;
        }
    }

    fn heap_up(&mut self, mut i: usize) {
        while i > 0 {
            let p = bin_parent(i);
            if self.heap_compare(p, i) {
                break;
            }
            self.heap_swap(p, i);
            i = p;
        }
    }

    // Small take/drop API to reuse nodes.
    fn node_take(
        &mut self, node_data: Node<TOrd, TData>,
    ) -> NodeHandle {
        let nhandle;
        if unlikely!(self.free == INVALID) {
            nhandle = self.node.len();
            self.node.push(node_data);
        } else {
            nhandle = self.free;
            let node = &mut self.node[nhandle];
            self.free = node.index;
            *node = node_data;
        }

        if cfg!(debug_assertions) {
            debug_assert!(self.contains_node_handle(&NodeHandle(nhandle)) == false);
        }

        NodeHandle(nhandle)
    }

    fn node_drop(
        &mut self, free_node: usize,
    ) -> TData {
        let node = &mut self.node[free_node];
        let user_data = node.user_data;
        node.index = self.free;
        self.free = free_node;
        user_data
    }

    // -------------------------------------------------------------------
    // Public API
    //
    pub fn insert(
        &mut self, value: TOrd, user_data: TData,
    ) -> NodeHandle {
        let tree_index = self.tree_index.len();

        let nhandle = self.node_take(Node {
            user_data: user_data,
            value: value,
            index: tree_index,
        });


        let index = self.tree_index.len();
        self.tree_index.push(nhandle.0);

        self.heap_up(index);

        // index in the self.nodes
        nhandle
    }

    pub fn pop_min(
        &mut self,
    ) -> Option<TData> {
        if self.tree_index.len() == 0 {
            return None;
        }

        let free_node = self.tree_index[0];

        if cfg!(debug_assertions) {
            debug_assert!(self.contains_node_handle(&NodeHandle(free_node)) == true);
        }

        let tree_index_len = self.tree_index.len() - 1;
        if tree_index_len != 0 {
            self.heap_swap(0, tree_index_len);
            self.tree_index.pop();
            self.heap_down(0);
        } else {
            self.tree_index.pop();
        }

        Some(self.node_drop(free_node))
    }

    pub fn pop_min_with_value(
        &mut self,
    ) -> Option<(TOrd, TData)> {
        // copied from pop_min
        if unlikely!(self.tree_index.len() == 0) {
            return None;
        }

        let free_node = self.tree_index[0];

        if cfg!(debug_assertions) {
            debug_assert!(self.contains_node_handle(&NodeHandle(free_node)) == true);
        }

        let tree_index_len = self.tree_index.len() - 1;
        if tree_index_len != 0 {
            self.heap_swap(0, tree_index_len);
            self.tree_index.pop();
            self.heap_down(0);
        } else {
            self.tree_index.pop();
        }
        // end copy from pop_min

        Some((self.node[free_node].value, self.node_drop(free_node)))
    }

    pub fn insert_or_update(
        &mut self, nhandle_p: &mut NodeHandle, value: TOrd, user_data: TData,
    ) {
        if nhandle_p.0 == INVALID {
            *nhandle_p = self.insert(value, user_data)
        } else {
            self.node_value_update_with_data(*nhandle_p, value, user_data);
        }
    }

    pub fn remove(&mut self, nhandle: NodeHandle) {
        debug_assert!(self.tree_index.len() != 0);
        debug_assert!(nhandle.0 < self.node.len());
        let mut i = self.node[nhandle.0].index;
        while i > 0 {
            let p = bin_parent(i);

            self.heap_swap(p, i);
            i = p;
        }
        self.pop_min();
    }

    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.tree_index.is_empty()
    }

    pub fn node_value_update(
        &mut self, nhandle: NodeHandle, value: TOrd,
    ) {
        debug_assert!(self.tree_index.len() != 0);
        debug_assert!(nhandle.0 < self.node.len());
        let (index, value_curr) = {
            let node = &mut self.node[nhandle.0];
            (node.index, node.value)
        };
        if value < value_curr {
            self.node[nhandle.0].value = value;
            self.heap_up(index);
        } else if value > value_curr {
            self.node[nhandle.0].value = value;
            self.heap_down(index);
        }
    }

    pub fn node_value_update_with_data(
        &mut self, nhandle: NodeHandle, value: TOrd, user_data: TData,
    ) {
        debug_assert!(self.tree_index.len() != 0);
        debug_assert!(nhandle.0 < self.node.len());
        self.node[nhandle.0].user_data = user_data;
        self.node_value_update(nhandle, value);
    }

    #[allow(dead_code)]
    pub fn node_value(
        &self, nhandle: NodeHandle,
    ) -> TOrd {
        self.node[nhandle.0].value
    }
    #[allow(dead_code)]
    pub fn node_data(
        &self, nhandle: NodeHandle,
    ) -> TData {
        self.node[nhandle.0].user_data
    }

    #[allow(dead_code)]
    pub fn new() -> MinHeap<TOrd, TData> {
        MinHeap {
            tree_index: vec![],
            node: vec![],
            free: INVALID,
        }
    }

    pub fn with_capacity(
        capacity: usize,
    ) -> MinHeap<TOrd, TData> {
        MinHeap {
            tree_index: Vec::with_capacity(capacity),
            node: Vec::with_capacity(capacity),
            free: INVALID,
        }
    }
}
