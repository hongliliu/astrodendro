class Leaf(object):

    ###########################################################################
    #   The following methods are used only by the Dendrogram class,          #
    #   for computing the dendrogram and should never be used manually:       #
    ###########################################################################
    def __init__(self, coord, f, idx=None):
        if not hasattr(f, '__iter__'): # Normal initialization - coord and f are for a single pixel:
            self.coords = [coord]
            self.f = [f]
            self.fmin, self.fmax = f, f
        else:
            self.coords = coord
            self.f = f
            self.fmin, self.fmax = min(f), max(f)
        self.idx = idx
        self.parent = None
        self._ancestor = None # Cached ancestor, if any
        self._level = None # Cached "level" property - see below

    def add_point(self, coord, f):
        "Add point to current leaf"
        self.coords.append(coord)
        self.f.append(f)
        self.fmin, self.fmax = min(f, self.fmin), max(f, self.fmax)

    def merge(self, leaf):
        self.coords.extend(leaf.coords)
        self.f.extend(leaf.f)
        self.fmin, self.fmax = min(leaf.fmin, self.fmin), max(leaf.fmax, self.fmax)

    def add_footprint(self, image, level):
        "Fill in a map which shows the depth of the tree"
        for c in self.coords:
            image[c] = level

    ###########################################################################
    #   The following methods can be used during OR after computation         #
    ###########################################################################
    
    @property
    def npix(self):
        return len(self.f)
    @property
    def npix_self(self):
        return len(self.f)
    @property
    def f_sum(self):
        return sum(self.f)
    @property
    def f_sum_self(self):
        return sum(self.f)
    @property
    def height(self):
        if self.parent == None:
            return self.fmax - self.fmin
        else:
            return self.fmax - self.parent.merge_level
    
    @property
    def ancestor(self):
        """
        Find the ancestor of this leaf/branch recursively. Always returns an
        object (may return self if the object has no parent).
        Results are partially cached to reduce recursion depth.
        The caching assumes that once an object has been given a parent, that
        parent will never change.
        
        This method should be equivalent to:
            if self.parent == None:
                return self
            else:
                return self.parent.ancestor
        """
        if self.parent == None:
            return self
        if not self._ancestor:
            self._ancestor = self.parent
        while self._ancestor.parent:
            # Use a loop rather than recursion to update
            # the cached ancestor, if needed:
            a = self._ancestor
            if a._ancestor:
                self._ancestor = a._ancestor # Update our cached value
            else:
                self._ancestor = a.parent
        return self._ancestor
    
    ###########################################################################
    #   The following methods are only reliable after the entire tree is      #
    #   computed. They should not be used in dendrogram.py                    #
    ###########################################################################

    @property
    def level(self):
        " Level: 0 for items in the trunk, 1 for their immediate children, etc"
        if self._level is None:
            if not self.parent:
                self._level = 0
            elif self.parent._level is not None:
                self._level = self.parent._level + 1
            else:
                #We could just use:
                #  self._level = self.parent.level + 1
                #But to avoid recursion, and keep things fast, we do it this way instead:
                obj = self.parent
                diff = 1
                while obj._level is None:
                    obj = obj.parent
                    diff += 1
                    # Note: we are counting on the dendrogram computation to ensure that
                    # the _level property of all items in the trunk has been set to zero
                self._level = obj._level + diff
                self.parent._level = self._level - 1
        return self._level

    @property
    def newick(self):
        " Newick representation of this Leaf " 
        return "%i:%.3f" % (self.idx, self.height)

    @property
    def peak(self):
        " Return (coordinates, flux) for the pixel with maximum value "
        self._peak = getattr(self, '_peak', (self.coords[self.f.index(self.fmax)], self.fmax)) 
        return self._peak 
    
    def get_peak_recursive(self):
        """
        Return the item, coordinate, and flux value of the pixel with highest
        intensity found within this leaf.
        This method exactly matches the more detailed method on Branch()
        """
        c,f = self.peak
        return self, c, f
    
    @property
    def eccentricity(self):
        """ Compute the max distance between this item and all other connected items """
        if type(self) is Leaf and self.parent == None:
            return 0
        max_distance = 0
        # First check all children of this item:
        if type(self) is Branch and self.items:
            max_distance = max([d.level for d in self.descendants]) - self.level
        # Now check other descendants of this item's ancestor
        ignore_item = self
        current_item = self.parent
        while current_item:
            peer_levels = [(max([d.level for d in item.descendants]) if type(item) is Branch else 1)
                           for item in current_item.items if item is not ignore_item]
            if peer_levels:
                distance = max(peer_levels) - current_item.level + (self.level - current_item.level)
                if distance > max_distance:
                    max_distance = distance
            ignore_item = current_item
            current_item = current_item.parent
        return max_distance

class Branch(Leaf):

    def __init__(self, items, coord, f, idx=None):
        self.merge_level = f # Record the exact flux level that triggered creation of this branch
        self.items = items
        for item in items:
            item.parent = self
        Leaf.__init__(self, coord, f, idx=idx)
        self._children_peak_result = None # cached result of get_peak_recursive for child with highest peak
        self._descendants = None # Cached value is not initially set

    ###########################################################################
    #   The following methods can be used during OR after computation         #
    ###########################################################################

    @property
    def npix(self):
        " Number of pixels among this branch and it's descendants "
        return len(self.f) + self.npix_children
    @property
    def npix_children(self):
        " Number of pixels among this branch's children "
        return sum([item.npix for item in self.items])
    @property
    def npix_self(self):
        " Number of pixels in this Branch, not counting children "
        return len(self.f)
    
    @property
    def f_sum(self):
        return sum(self.f) + self.f_sum_children
    @property
    def f_sum_children(self):
        return sum([item.f_sum for item in self.items])
    @property
    def f_sum_self(self):
        return sum(self.f) 
    
    def add_footprint(self, image, level, recursive=True):
        if recursive:
            for item in self.items:
                item.add_footprint(image, level + 1)
        Leaf.add_footprint(self, image, level)
    
    
    ###########################################################################
    #   The following methods are only reliable after the entire tree is      #
    #   computed. They should not be used in dendrogram.py                    #
    ###########################################################################

    @property
    def newick(self):
        newick_items = [item.newick for item in self.items]
        return "(%s)%s:%.3f" % (','.join(newick_items), self.idx, self.height)
    
    @property
    def descendants(self):
        "Returns a flattened list of all child leaves and branches. Non-recursive."
        if self._descendants is None:
            self._descendants = []
            items_to_add = [self] # More specifically, items whose children we will need to add to the list
            while True:
                children = []
                map(children.extend, [i.items for i in items_to_add])
                self._descendants.extend(children)
                # Then proceed, essentially recursing through child Branch items:
                items_to_add = [b for b in children if type(b) is Branch]
                if not items_to_add:
                    break
        return self._descendants
    
    def get_peak_recursive(self):
        """
        Return the item, coordinate, and flux value of the pixel with highest
        intensity among this branch and its children.
        """
        if self._children_peak_result: # If we have cached the result of peak flux among children:
            ci, cc, cf = self._children_peak_result
        else:
            # get Child Item, Child peak Coords, Child Flux for child with highest max. flux
            ci, cc, cf = self.items[0].get_peak_recursive()
            for child in self.items[1:]:
                ci_try, cc_try, cf_try = child.get_peak_recursive() 
                if cf_try > cf:
                    ci, cc, cf = ci_try, cc_try, cf_try
            self._children_peak_result = ci, cc, cf
        if cf >= self.fmax: # A child has a higher flux than this branch does:
            return ci, cc, cf
        else:
            # This branch directly contains a pixel with higher flux than any
            # children have. This happens sometimes if min_delta and/or
            # min_npix are nonzero
            return self, self.coords[self.f.index(self.fmax)], self.fmax