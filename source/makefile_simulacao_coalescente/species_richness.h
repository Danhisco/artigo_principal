/************************************************************
  DOCUMENTATION
 ************************************************************
 APPENDIX TO
 "A COALESCENCE APPROACH TO SPATIAL NEUTRAL ECOLOGY"
 James Rosindell, Yan Wong, Rampal Etienne
 (release version 1.0)
License Agreement:
This software is provided free of charge.
You may use and distribute this code freely.
By using this code you agree to correctly
cite the manuscript
"A COALESCENCE APPROACH TO SPATIAL NEUTRAL ECOLOGY"
James Rosindell, Yan Wong, Rampal Etienne

in any publications that utilize this code.
You use this code at your own risk.
*/
# include <vector>
/************************************************************
  ARRAY DATA STRUCTURE OBJECT
 ************************************************************/
// Template class array1D
template<class X>
class array1D {
    public:
        array1D();
        ~array1D();
        void SetSize(int n);

        int size();
        X& operator[](int index);
    private:
        int entries; // stores the number of entries
        X* data; // stores the entries themselves
};
/************************************************************
  TREE OBJECT
 ************************************************************/
class tree
// this object represents the output coalescence tree itself
// and has all the useful functions - everything above this point is just a tool
// that is required for this object
// the end user should initialise one instance of "tree" and
// use that object to do all their calculations
{
    public:
        tree();
        //PRODUCING A COALESCENCE TREE IN CLASS TREE
        // this is the part of the code that actually creates the colescence tree
        // input variables are described in detail below
        // area1 = width of survey area (in trees)
        // area2 = length of survey area (int trees)
        // minspec = smallest speciation rate required
        // dispersal = dispersal distance (dispersal kernel width)
        // tol = tolerance in results (0.01 means 1% error - 0.5% either side)
        // typeflag deals with the type of dispersal in the model (true means normal)
        //void maketree(int area1 , int area2 , double minspec , int dispersal , double tol, bool normflag);

        // this returns an interval within which the true mean richness must lie
        array1D<double> get_richnessint(double spec);
        // this returns the midpoint between the maximum and minimum richness estimates
        double get_richness(double spec);
        class treenode
            // this is a data storage object for the coalescence tree
            // an array of these objects is required to store the whole tree
            // we include a few simple functions that prove useful later
        {
            public:
                treenode();
                void setup(bool root_);
                void set_root(bool root_);
                void set_parent(unsigned int x);
                void set_steps(unsigned int x);
                bool get_root();
                unsigned int get_parent();
                int get_steps();
                // increments the steps by one
                void inc_steps();
            // violation of privacy!? :p
            //private:
                bool root;
                // is this node at the end of the tree (true)
                // or just here to mark a coalescense (false)
                unsigned int parent;
                // this stores the parent of the individual
                // 0 means there is no parent - we are at the end of the tree
                // (as far as has been calculated)
                unsigned int steps;
                // the number of generations (chances of speciation)
                // between this individual and its parent in the tree
        };
    // converts array of labels and steps to tree structure
    //void convert_from_labelarray(double minspec, int size_pilha, int *v_label, int *coalescence_events, int *steps);
    // violation of privacy!? :p
    //private:
        array1D<treenode> data;
        // stores the coalescence tree itself
        unsigned int enddata;
        // marks the end of the data array
        // this is so that size can be expanded and contracted easily up to a maximum size
        double minspecsetup;
        // when producing a coalescence tree, we do so with a minimal speciation rate
        // in mind (to save a lot of computational expense in growing the full tree)
        // when true some messages will be printed out for debub purposes

        unsigned int get_greatest_parent(unsigned int x);
        // gets the greatest parent of node x

        std::vector< std::vector<int> > children;
        // this vector stores the extra information needed for Newick outputs
        // it is left clear as standard and only used when needed
};

