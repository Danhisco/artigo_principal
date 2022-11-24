# include <stdio.h>
# include <fstream>
# include <vector>
# include <iostream>
# include <string>
# include <math.h>

using namespace std;
#include "species_richness.h"
/************************************************************
  ARRAY DATA STRUCTURE OBJECT
 ************************************************************/
// Template class array1D
    template<class X>
array1D<X>::array1D():data(0)
{
    SetSize(0);
}
// destructor
    template<class X>
array1D<X>::~array1D()
{

    delete[] data;
}
// setter for size
    template<class X>
void array1D<X>::SetSize(int n)
{
    if(data)
    {
        delete[] data;
    }
    if (n>0)
    {
        data=new X[n];
    }
    else
    {
        data=0;
    }
    entries=n;
}
// getter for size
    template<class X>
int array1D<X>::size()
{
    return entries;
}
// [] operator
    template<class X>
X& array1D<X>::operator[](int index)
{
    if (index < 0)
    {
        index = -1 * index;
        index = index % entries;

        index = entries - index;
    }
    index = index % entries;
    return data[index];
}

/************************************************************
  TREE OBJECT
 ************************************************************/
tree::treenode::treenode()
{
}
void tree::treenode::setup(bool root_)
    // sets up variables with initial conditions
{
    root = root_;
    parent = 0;
    steps = 0;
}

// standard setters
void tree::treenode::set_root(bool root_)
{
    root = root_;
}
void tree::treenode::set_parent(unsigned int x)
{
    parent = x;
}
void tree::treenode::set_steps(unsigned int x)
{
    steps = x;
}
// standard getters
bool tree::treenode::get_root()
{
    return root;
}
unsigned int tree::treenode::get_parent()
{
    return parent;
}
unsigned int tree::get_greatest_parent(unsigned int x)
{
    while (data[x].get_parent())
        x = data[x].get_parent();
    return x;
}
int tree::treenode::get_steps()
{
    return (steps);
}
// increments the steps by one
void tree::treenode::inc_steps()
{
    steps++;
}
/************************************************************
  CLASS TREE - THE MAIN CLASS
 ************************************************************/
tree::tree()
    // initialisor
{
    minspecsetup = 2;
    children.clear();
}
// TODO Is this even necessary?
//explicitly instantiate array1D for treenodes and for doubles
template class array1D<tree::treenode>;
template class array1D<double>;

/************************************************************
  RICHNESS CALCULATION METHODS
 ************************************************************/
array1D<double> tree::get_richnessint(double spec)
    // this returns an interval within which the true mean ricness must lie
{
    array1D<double> result;
    result.SetSize(2);
    result[0] = 0.0;
    result[1] = 0.0;
    if (minspecsetup <= spec)
        // check that the tree was calculated for a small enough speciation rate
        // it is possible to override this check by commenting out if required
    {
        // probarray stores, for each node, a probability that is required in the calculation
        array1D<double> probarray;
        probarray.SetSize(enddata+1);
        for (int i = 1 ; i <= enddata ; i ++)
            // loop over all nodes and initialise the relating element in probarray
        {
            if (data[i].get_root() == false)
                // a value of -2.0 indicates an internal node
                // that has thus far not been pruned at all
            {
                probarray[i] = -2.0;
            }
            else
                // a value of 1.0 is assigned to every free branch
                // this is because it is certain that the lineages have not encountered speciaiton

                // when they are at the very end and no pruning has so far taken place
            {
                probarray[i] = 1.0;
            }
        }
        bool loop = true;
        while (loop)
            // continue looping until we have completed our calculations
        {
            loop = false;
            for (int i = 1 ; i <= enddata ; i++)
            {
                // check to see if that part of the array is complete
                if (probarray[i] >= 0.0)
                {
                    // it is complete so use it to complete the info on in its parents
                    int theparent = data[i].get_parent();
                    if ((probarray[theparent] < 0.0) && (theparent != 0))
                    {
                        // only do anything if the parents have not already been completed
                        loop = true;
                        // be sure to go round again if the parents have not been completed
                        if (probarray[theparent] <= -1.5)
                        {
                            // parent not at all complete
                            double temprob = (pow(1-spec,double(data[i].get_steps())));
                            result[0] += probarray[i]*(1-temprob);
                            // we store probabilities as negative in a node if they
                            // refer only to one of the two branches of the node
                            // we then wait until we have both branches of the node
                            // before continuing to complete the full calculation
                            probarray[theparent]= -1.0*probarray[i]*temprob;
                        }
                        else
                        {
                            // parent partailly complete
                            double temprob = (pow(1-spec,double(data[i].get_steps())));

                            // update Smin
                            result[0] += probarray[i]*(1-temprob);
                            // update the probability array
                            temprob = temprob*probarray[i];
                            double temprob2 = probarray[theparent]*-1.0;
                            probarray[theparent]=(temprob+temprob2-temprob*temprob2);
                        }
                    }
                }
                else
                {
                    // be sure to repeat the loop unless all calculations are fully compelted
                    loop = true;
                }
            }
        }
        for (int i = 1 ; i <= enddata ; i++)
        {
            // here we are dealing with all the last branches after prooning all nodes
            if (data[i].get_parent() == 0)
            {
                result[0] += probarray[i]*(1-pow(1-spec,double(data[i].get_steps())));
                result[1] += probarray[i]*pow(1-spec,double(data[i].get_steps()));
            }
        }
        // return the result
        result[1] += result[0];
        return result;
    }
    else
    {
        // return a default error result to indicate
        // the tree was unsutable for the requested calculations
        result[0] = -1.0;
        result[1] = -1.0;
        return result;
    }

}

double tree::get_richness(double spec)
    // this returns the midpoint between the maximum and minimum richness estimates
{
    array1D<double> richnessresult = get_richnessint(spec);
    return (richnessresult[0]+richnessresult[1])/2.0;
}


/* 
 * converts array of labels and steps to tree structure
 */
/*void tree::convert_from_labelarray(double minspec, int size_pilha, int *v_label, int *coalescence_events, int *steps)
{
    int origin, chosen;
    minspecsetup = minspec;
    // data - this will store the coalescence tree itself
    // there can only be a maximum of twice as many nodes as there are
    // initially free branches so we can set the size of our data object
    data.SetSize(2*size_pilha);
    enddata = 0; // 0 is reserved as null
    // this means that data is empty and that data[1] is where the first
    // piece of data will be recorded
    for (int i = 0 ; i < size_pilha; i++)
    {
        // looping over the entire survey area
        // setup the data variable
        enddata ++;
        data[enddata].setup(true); // generation number = 1
    }

    for (int i = 0; i < size_pilha; i++){
        origin = coalescence_events[i];
        chosen = v_label[i];
        while (chosen < 0){
            current = get_greatest_parent(chosen-1);
            enddata++;
            data[enddata].setup(false);
            chosen = -chosen;
            data[origin-1].set_parent(enddata);
            data[current].set_parent(enddata);
            data[origin-1].set_steps(steps[origin]);
            data[current].set_steps(steps[chosen]);
            origin = chosen;
            chosen = v_label[chosen];
            if (chosen < 0)
                data[enddata].setup(true);
        }
    }
}*/
