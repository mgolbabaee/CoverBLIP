#include <vector>
#include <algorithm>

#include<string.h>
#include<math.h>
#include<list>
#include<stdlib.h>
//#define NDEBUG
#include<assert.h>
#include <limits.h>
#include <values.h>
#include <stdint.h>
#include <iostream>

#include<omp.h>

#include "clspoint.cc"

struct node {
  point p;
  double max_dist;  // The maximum distance to any grandchild.
  double parent_dist; // The distance to the parent.
  node* children;
  unsigned short int num_children; // The number of children.
  short int scale; // Essentially, an upper bound on the distance to any child.
};

struct ds_node{
  v_array<double> dist;
  point p;
};

struct d_node {
  double dist;
  const node *n;
};

#define SWAP(a, b)				\
  do						\
    {						\
      d_node tmp = * a;				\
      * a = * b;				\
      * b = tmp;				\
    } while (0)

#define  TR  covertree
class TR
{
public:
    TR (); // constructor
    ~TR (); // DE-structor
    
    // initializing
    void Get_Dict_From_Matrix(const mxArray* mxData);
    void Get_Epsilon_Aprox_From_Matrix(const mxArray* mxData);
    void Get_Stop_Scale_From_Matrix(const mxArray* mxData);
    int Get_Num_Atom_From_Obj();
    int Get_Max_Tree_Level_From_Obj();

    // create/destroy
    void Destroy_Tree();
    void Build_Tree();
    void Print_Tree();
    void print_tree_structure(mxArray* mxOutput_tree_structure);
    void eNN_loop_current(const mxArray* mxInput_uppbound, const mxArray* mxInput_query, mxArray* mxOutput_label, mxArray* mxOutput_atom, mxArray* mxOutput_dmins, mxArray* mxOutput_counter);
    void eNN_bsearch(const mxArray* mxInput_batchsize, const mxArray* mxInput_query, mxArray* mxOutput_label, mxArray* mxOutput_atom, mxArray* mxOutput_dmins, mxArray* mxOutput_counter, mxArray* mxOutput_index);
    void eNN_stop_scale(const mxArray* mxInput_query, mxArray* mxOutput_label, mxArray* mxOutput_atom, mxArray* mxOutput_dmins, mxArray* mxOutput_counter);

private:
   v_array<point> dict;
   node root;
   double internal_epsilon_approx;
   int counter;
   int dict_num_atom;
   int max_tree_level;
   int stop_scale;
   clspoint*  Point_Obj;

    void print_tree(int depth, node &top_node);
    void free_tree(node &top_node);
    void tree_structure(node &top_node, double* tree_table);
 
    //information gathering
    int height_dist(const node top_node,v_array<int> &heights);
    void breadth_dist(const node top_node,v_array<int> &breadths);
    void depth_dist(int top_scale, const node top_node,v_array<int> &depths);   
    
    float base = 1.3;
    float il2 = 1./log(base);
    inline double dist_of_scale(int s);
    inline int get_scale(double d);
    int min(int f1, int f2);
    double max(double f1, double f2);
    node new_node(const point &p);
    node new_leaf(const point &p);
    double max_set(v_array<ds_node> &v);
    
    void print_space(int s);
    void split(v_array<ds_node>& point_set, v_array<ds_node>& far_set, int max_scale);
    void dist_split(v_array<ds_node>& point_set, v_array<ds_node>& new_point_set, point new_point, int max_scale);
    node batch_insert(const point& p, int max_scale, int top_scale,v_array<ds_node>& point_set, v_array<ds_node>& consumed_set,v_array<v_array<ds_node> >& stack);
    node batch_create(v_array<point> points);   
    void add_height(int d, v_array<int> &heights); 

    inline double compare(const d_node *p1, const d_node* p2);
    void halfsort (v_array<d_node > cover_set);
    v_array<v_array<d_node> > get_cover_sets(v_array<v_array<v_array<d_node> > > &spare_cover_sets);
    inline bool shell(double parent_query_dist, double child_parent_dist, double upper_bound);
    
    void update(double *k_upper_bound, double upper_bound);
    double *alloc_upper();
    void setter(double* begin, double max);
    
    inline void copy_zero_set(node* query_chi, double* new_upper_bound, v_array<d_node> &zero_set, v_array<d_node> &new_zero_set);
    inline void copy_cover_sets(node* query_chi, double* new_upper_bound, v_array<v_array<d_node> > &cover_sets,v_array<v_array<d_node> > &new_cover_sets, int current_scale, int max_scale);    
    void print_query(const node *top_node);
    void print_cover_sets(v_array<v_array<d_node> > &cover_sets, v_array<d_node> &zero_set, int current_scale, int max_scale);
    inline void descend(const node* query, double* upper_bound, int current_scale, int &max_scale, v_array<v_array<d_node> > &cover_sets, v_array<d_node> &zero_set);
    inline void descend_approx(const node* query, double* upper_bound, int current_scale, int &max_scale, v_array<v_array<d_node> > &cover_sets, v_array<d_node> &zero_set);
    void brute_nearest(const node* query,v_array<d_node> zero_set, double* upper_bound, v_array<v_array<point> > &results, v_array<double> &dist_mins, v_array<v_array<d_node> > &spare_zero_sets);
    void internal_batch_nearest_neighbor(const node *query, v_array<v_array<d_node> > &cover_sets,v_array<d_node> &zero_set,int current_scale,int max_scale,double* upper_bound,v_array<v_array<point> > &results, v_array<double> &dist_mins, v_array<v_array<v_array<d_node> > > &spare_cover_sets,v_array<v_array<d_node> > &spare_zero_sets);
    void batch_nearest_neighbor(const node &top_node, const node &query, v_array<v_array<point> > &results, v_array<double> &dist_mins);
    void batch_nearest_neighbor_current(const node &top_node, const node &query, v_array<v_array<point> > &results, v_array<double> &dist_mins, double input_uppbound);
    void epsilon_approx_nearest_neighbor(const node &top_node, const node &query, v_array<v_array<point> > &results, double epsilon_approx);

};

/***************************************************************************************/
/* constructor */
TR::TR ()
{
    // init
    internal_epsilon_approx = 0.;
    stop_scale = 0;
    counter = 0;
    max_tree_level = 150;
    //dmin = MAXFLOAT;
    Point_Obj = new clspoint;
}
/***************************************************************************************/


/***************************************************************************************/
/* destructor */
TR::~TR ()
{
    // delete the tree
    Destroy_Tree();
}
/***************************************************************************************/

/***************************************************************************************/
/* destroy! */
void TR::Destroy_Tree()
{
    Point_Obj->freeup(dict);
    free_tree(root);
}
/***************************************************************************************/

void TR::Get_Dict_From_Matrix(const mxArray* mxData)
{
    int Pm, Pn;
    Pm = mxGetM(mxData);
    Pn = mxGetN(mxData);
    float *Pr, *Pi;
    Pr = (float *)mxGetPr(mxData);
    if(mxIsComplex(mxData))
        Pi = (float *)mxGetPi(mxData); 
    else
        Pi = 0; 
    dict = Point_Obj->parse_points(Pr,Pi, Pm, Pn);
    dict_num_atom = Pm;
 }

int TR::Get_Num_Atom_From_Obj()
{
    return dict_num_atom;
}

int TR::Get_Max_Tree_Level_From_Obj()
{
    return max_tree_level;
}

void TR::Get_Epsilon_Aprox_From_Matrix(const mxArray* mxData)
{
    double *ep;
    ep = mxGetPr(mxData);
    internal_epsilon_approx = ep[0]; 
}

void TR::Get_Stop_Scale_From_Matrix(const mxArray* mxData)
{
    double *ep;
    ep = mxGetPr(mxData);
    stop_scale = ep[0]; 
}

void TR::Build_Tree()
{
    root = batch_create(dict);  
  //  mexPrintf("After dictionary generation: Counter =  %d!\n", counter);
    counter = 0;
}

void TR::print_tree_structure(mxArray* mxOutput_tree_structure)
{
    double *tree_table;
    tree_table = mxGetPr(mxOutput_tree_structure);
    tree_structure(root, tree_table);
    //print_tree(5, root);
}

void TR::tree_structure(node &top_node, double* tree_table)
{
    int level;
    if(top_node.children>0)
    {
	    node *chi = top_node.children;
	    node *child_end = top_node.children + top_node.num_children;
        int i = 0;
	    for (chi; chi != child_end; chi++)
        {
            level =chi->scale;
            tree_table[(chi->p->label-1)*max_tree_level + level-1] = (double) top_node.p->label;
            tree_structure(top_node.children[i], tree_table);
            i =i+1;
        }
    }
}

void TR::free_tree(node &top_node)
{
    if(top_node.children>0)
    {
        for(int i = 0; i < top_node.num_children; i++)
            free_tree(top_node.children[i]);
        free(top_node.children);
    }
    return;
}

void TR::eNN_loop_current(const mxArray* mxInput_uppbound,const mxArray* mxInput_query, mxArray* mxOutput_label, mxArray* mxOutput_atom, mxArray* mxOutput_dmins, mxArray* mxOutput_counter)
{
    int Qm, Qn;
    Qm = mxGetM(mxInput_query);
    Qn = mxGetN(mxInput_query);
    float *Qr, *Qi;
    Qr = (float *)mxGetPr(mxInput_query);
    if(mxIsComplex(mxInput_query))
        Qi = (float *)mxGetPi(mxInput_query);
    else
        Qi = 0;
    double *U;
    U = mxGetPr(mxInput_uppbound);

    double *Rl, *Ra_r, *Ra_i, *Rd, *Rc;
    Rl = mxGetPr(mxOutput_label);
    Ra_r = mxGetPr(mxOutput_atom);
    Ra_i = mxGetPi(mxOutput_atom);
    Rd = mxGetPr(mxOutput_dmins);
    Rc = mxGetPr(mxOutput_counter);

    v_array<point> set_of_queries = Point_Obj->parse_points(Qr,Qi, Qm, Qn); 
    v_array<v_array<point>> res;    
    v_array<double> dmins;
    v_array<point> temp;
    #pragma omp parallel for private(res, dmins) num_threads(0)
        for(int i=0; i<Qm; i++)
        {
            batch_nearest_neighbor_current(root,new_leaf(set_of_queries[i]),res, dmins, U[i]);
            Rl[i] = 0;
            if(res[0].index > 1)
            {
                Rl[i] = res[0][1][0].label;
                Point_Obj->output(Ra_r+i*Qn, Ra_i+i*Qn,res[0][1]);
            }
            for(int j=0; j<res[0].index; j++)
                temp = pop(res);
            Rd[i] = dmins[0]; 
            dmins.index = 0;  //
            res.index = 0;
        }  
    Point_Obj->freeup(set_of_queries);
    Point_Obj->freeup2(res);
    Rc[0] = counter;
    counter = 0;
}

void TR::eNN_bsearch(const mxArray* mxInput_batchsize,
const mxArray* mxInput_query, 
mxArray* mxOutput_label, mxArray* mxOutput_atom, 
mxArray* mxOutput_dmins, mxArray* mxOutput_counter, 
mxArray* mxOutput_index)
{
    int Qm, Qn;
    Qm = mxGetM(mxInput_query);
    Qn = mxGetN(mxInput_query);
    float *Qr, *Qi;
    Qr = (float *)mxGetPr(mxInput_query);
    if(mxIsComplex(mxInput_query))
        Qi = (float *)mxGetPi(mxInput_query);
    else
        Qi = 0;
    double *B;
    B = mxGetPr(mxInput_batchsize);

    double *Rl, *Ra_r, *Ra_i, *Rd, *Rc, *Ri;
    Rl = mxGetPr(mxOutput_label);
    Ra_r = mxGetPr(mxOutput_atom);
    Ra_i = mxGetPi(mxOutput_atom);
    Rd = mxGetPr(mxOutput_dmins);
    Rc = mxGetPr(mxOutput_counter);
    Ri = mxGetPr(mxOutput_index);

    v_array<point> set_of_queries; // = Point_Obj->parse_points(Qr,Qi, Qm, Qn); 
    node query_root; // = batch_create(set_of_queries);  
    v_array<v_array<point>> res;
    v_array<point> temp;
    v_array<double> dmins;
    int chop_size = B[0];
    int chop_num = Qm/chop_size;
    if(Qm%chop_size>0)
    {
         mexPrintf("For the moment, the number of query atoms divided by the batch size need to be an integer!");
         return;
     }
    //printf("maximum threads: %d\n", omp_get_max_threads());
    #pragma omp parallel for private(set_of_queries, query_root, res, dmins) num_threads(0)
        for(int i=0; i<chop_num; i++)
        {
            set_of_queries = Point_Obj->parse_points_chop(Qr+i*chop_size, Qi+i*chop_size, chop_size, Qm, Qn);
            query_root = batch_create(set_of_queries);
            batch_nearest_neighbor_current(root,query_root,res, dmins, MAXFLOAT);
           //printf("thread %d out of %d threads\n", omp_get_thread_num(), omp_get_num_threads());
           for(int j=0;j<chop_size;j++)
           {
            Rl[j+i*chop_size] = 0;
            Rd[j+i*chop_size] = dmins[j];
                if(res[j].index > 1)
                {
                    Rl[j+i*chop_size] = res[j][1][0].label;
                    Ri[j+i*chop_size] = res[j][0][0].label;
                    Point_Obj->output(Ra_r+i*chop_size+j*Qn, Ra_i+i*chop_size+j*Qn,res[j][1]);
                }
            }
           for(int k=0; k<res.index; k++)
                  temp = pop(res);
           res.index = 0;
           set_of_queries.index = 0;
           free_tree(query_root);
           dmins.index = 0;
        }   
    Point_Obj->freeup(set_of_queries);
    Rc[0] = counter;
    counter = 0;
}

void TR::eNN_stop_scale(const mxArray* mxInput_query, mxArray* mxOutput_label, mxArray* mxOutput_atom, mxArray* mxOutput_dmins, mxArray* mxOutput_counter)
{
    int Qm, Qn;
    Qm = mxGetM(mxInput_query);
    Qn = mxGetN(mxInput_query);
    float *Qr, *Qi;
    Qr = (float *)mxGetPr(mxInput_query);
    if(mxIsComplex(mxInput_query))
        Qi = (float *)mxGetPi(mxInput_query);
    else
        Qi = 0;

    double *Rl, *Ra_r, *Ra_i, *Rd, *Rc;
    Rl = mxGetPr(mxOutput_label);
    Ra_r = mxGetPr(mxOutput_atom);
    Ra_i = mxGetPi(mxOutput_atom);
    Rd = mxGetPr(mxOutput_dmins);
    Rc = mxGetPr(mxOutput_counter);

    v_array<point> set_of_queries = Point_Obj->parse_points(Qr,Qi, Qm, Qn); 
    v_array<v_array<point>> res;    
    v_array<double> dmins;
    v_array<point> temp;
    #pragma omp parallel for private(res, dmins) num_threads(0)
        for(int i=0; i<Qm; i++)
        {
            batch_nearest_neighbor_current(root,new_leaf(set_of_queries[i]),res, dmins, MAXFLOAT);
            Rl[i] = 0;
            if(res[0].index > 1)
            {
                Rl[i] = res[0][1][0].label;
                Point_Obj->output(Ra_r+i*Qn, Ra_i+i*Qn,res[0][1]);
            }
            for(int j=0; j<res[0].index; j++)
                temp = pop(res);
            Rd[i] = dmins[0]; 
            dmins.index = 0;  //
            res.index = 0;
        }  
    Point_Obj->freeup(set_of_queries);
    Point_Obj->freeup2(res);
    Rc[0] = counter;
    counter = 0;
}

inline double TR::dist_of_scale (int s)
{
  return std::pow(base, s);
}

inline int TR::get_scale(double d)
{
  return (int) ceilf(il2 * log(d));
}
#include<numeric>

int TR::min(int f1, int f2)
{
  if ( f1 <= f2 )
    return f1;
  else 
    return f2;
}

double TR::max(double f1, double f2)
{
  if ( f1 <= f2 )
    return f2;
  else 
    return f1;
}

node TR::new_node(const point &p)
{
  node new_node;
  new_node.p = p;
  return new_node;
}

node TR::new_leaf(const point &p)
{
  node new_leaf = {p,0.,0.,NULL,0,150};
  return new_leaf;
}

double TR::max_set(v_array<ds_node> &v)
{
  double max = 0.;
  for (int i = 0; i < v.index; i++)
    if ( max < v[i].dist.last()) 
      max = v[i].dist.last();
  return max;
}

void TR::print_space(int s)
{
  for (int i = 0; i < s; i++)
    mexPrintf(" ");
}

void TR::print_tree(int depth, node &top_node)
{
  print_space(depth);
  Point_Obj->print(top_node.p);
  if ( top_node.num_children > 0 ) {
    print_space(depth); mexPrintf("scale = %i\n",top_node.scale);
    print_space(depth); mexPrintf("max_dist = %f\n",top_node.max_dist);
    print_space(depth); mexPrintf("num children = %i\n",top_node.num_children);
    for (int i = 0; i < top_node.num_children;i++)
      print_tree(depth+1, top_node.children[i]);
  }
}

void TR::split(v_array<ds_node>& point_set, v_array<ds_node>& far_set, int max_scale)
{
  unsigned int new_index = 0;
  double fmax = dist_of_scale(max_scale);
  for (int i = 0; i < point_set.index; i++){
    if (point_set[i].dist.last() <= fmax) {
      point_set[new_index++] = point_set[i];
    }
    else
      push(far_set,point_set[i]);
  }
  point_set.index=new_index;  
}

void TR::dist_split(v_array<ds_node>& point_set, 
		v_array<ds_node>& new_point_set, 
		point new_point, 
		int max_scale)
{
  unsigned int new_index = 0;
  double fmax = dist_of_scale(max_scale);
  for(int i = 0; i < point_set.index; i++) 
    {
      double new_d;
      new_d = Point_Obj->distance(new_point, point_set[i].p, fmax, &counter);
      if (new_d <= fmax ) {
	push(point_set[i].dist, new_d);
	push(new_point_set,point_set[i]);
      }
      else
	point_set[new_index++] = point_set[i];
    }
  point_set.index = new_index;
}

/*
   max_scale is the maximum scale of the node we might create here.
   point_set contains points which are 2*max_scale or less away.
*/


node TR::batch_insert(const point& p, 
		  int max_scale, 
		  int top_scale,
		  v_array<ds_node>& point_set, 
		  v_array<ds_node>& consumed_set,
		  v_array<v_array<ds_node> >& stack)
{
  if (point_set.index == 0) 
    return new_leaf(p);
  else {
    double max_dist = max_set(point_set); //O(|point_set|)
    int next_scale = min (max_scale - 1, get_scale(max_dist));
    if (next_scale < -150) // We have points with distance 0.
      {
	v_array<node> children;
	push(children,new_leaf(p));
	while (point_set.index > 0)
	  {
	    push(children,new_leaf(point_set.last().p));
	    push(consumed_set,point_set.last());
	    point_set.decr();
	  }
	node n = new_node(p);
	n.scale = 150; // A magic number meant to be larger than all scales.  
	n.max_dist = 0;
	alloc(children,children.index);
	n.num_children = children.index;
	n.children = children.elements;
	return n;
      }
    else
      {
	v_array<ds_node> far = pop(stack);
	split(point_set,far,max_scale); //O(|point_set|)
	
	node child = batch_insert(p, next_scale, top_scale, point_set, consumed_set, stack);
	
	if (point_set.index == 0)
	  {
	    push(stack,point_set);
	    point_set=far;
	    return child;
	  }
	else {
	  node n = new_node(p);
	  v_array<node> children;
	  push(children, child);
	  v_array<ds_node > new_point_set = pop(stack);
	  v_array<ds_node > new_consumed_set = pop(stack);
	  while (point_set.index != 0) { //O(|point_set| * num_children)
	    point new_point = point_set.last().p;
	    double new_dist = point_set.last().dist.last();
	    push(consumed_set, point_set.last());
	    point_set.decr();
	    
	    dist_split(point_set, new_point_set, new_point, max_scale); //O(|point_saet|)
	    dist_split(far,new_point_set,new_point,max_scale); //O(|far|)
	    
	    node new_child = 
	      batch_insert(new_point, next_scale, top_scale, new_point_set, new_consumed_set, stack);
	    new_child.parent_dist = new_dist;

	    push(children, new_child);
	    
	    double fmax = dist_of_scale(max_scale);
	    for(int i = 0; i< new_point_set.index; i++) //O(|new_point_set|)
	      {
		new_point_set[i].dist.decr();
		if (new_point_set[i].dist.last() <= fmax)
		  push(point_set, new_point_set[i]);
		else
		  push(far, new_point_set[i]);
	      }
	    for(int i = 0; i< new_consumed_set.index; i++) //O(|new_point_set|)
	      {
		new_consumed_set[i].dist.decr();
		push(consumed_set, new_consumed_set[i]);
	      }
	    new_point_set.index = 0;
	    new_consumed_set.index = 0;
	  }
	  push(stack,new_point_set);
	  push(stack,new_consumed_set);
	  push(stack,point_set);
	  point_set=far;
	  n.scale = top_scale - max_scale;   
	  n.max_dist = max_set(consumed_set);
	  alloc(children,children.index);
	  n.num_children = children.index;
	  n.children = children.elements;
	  return n;
	}
      }
  }
}
  
node TR::batch_create(v_array<point> points) 
{
  assert(points.index > 0);
  v_array<ds_node > point_set;
  v_array<v_array<ds_node > > stack;

  for (int i = 1; i < points.index; i++) {
    ds_node temp;
    push(temp.dist, Point_Obj->distance(points[0], points[i], MAXFLOAT, &counter)); 
    temp.p = points[i];
    push(point_set,temp);
  }
  v_array<ds_node> consumed_set;
  
  double max_dist = max_set(point_set);

  node top = batch_insert (points[0], 
			   get_scale(max_dist),
			   get_scale(max_dist),
			    point_set, 
			    consumed_set,
			    stack);
  for (int i = 0; i<consumed_set.index;i++)
    free(consumed_set[i].dist.elements);
  free(consumed_set.elements);
  for (int i = 0; i<stack.index;i++)
    free(stack[i].elements);
  free(stack.elements);
  free(point_set.elements);
  return top;
}

void TR::add_height(int d, v_array<int> &heights)
{
  if (heights.index <= d)
    for(;heights.index <= d;)
      push(heights,0);
  heights[d] = heights[d] + 1;
}

int TR::height_dist(const node top_node,v_array<int> &heights)
{
  if (top_node.num_children == 0)
    {
      add_height(0,heights);
      return 0;
    }
  else
    {
      int max_v=0;
      for (int i = 0; i<top_node.num_children ;i++)
	{
	  int d = height_dist(top_node.children[i], heights);
	  if (d > max_v)
	    max_v = d;	
	}
      add_height(1 + max_v, heights);
      return (1 + max_v);
    }
}

void TR::depth_dist(int top_scale, const node top_node,v_array<int> &depths)
{
  if (top_node.num_children > 0)
      for (int i = 0; i<top_node.num_children ;i++)
	{
	  add_height(top_node.scale, depths);
	  depth_dist(top_scale, top_node.children[i], depths);
	}
}

void TR::breadth_dist(const node top_node,v_array<int> &breadths)
{
  if (top_node.num_children == 0)
    add_height(0,breadths);
  else
    {
      for (int i = 0; i<top_node.num_children ;i++)
	breadth_dist(top_node.children[i], breadths);
      add_height(top_node.num_children, breadths);
    }
}
/***************************************************************************************/
/* print the tree info to the MATLAB display as ASCII text. */
void TR::Print_Tree()
{
    mexPrintf("Start from the ROOT Level:\n\n");
    
    print_tree(5,root);

 //   print_count = 0;  // initialize (for error checking)
//    Print_Tree(root);

 /*   const UINT_type Num_Points = (UINT_type) Points.Get_Num_Rows();
  //  const UINT_type Num_Cols = (UINT_type) Points.Get_Num_Cols();
   // for(int i=0; i<Num_Points; i++)
    {
        //for(int j =0; j<Num_Cols; j++)
        { const pt_type* test;
          test = Points.Get_Data_Col_Ptr(1);
          mexPrintf("%f ", test[i]);}
        mexPrintf("\n");
    }*/
/*    if (print_count!=Num_Points)
        {
        mexPrintf("Number of points printed: %d.\n",print_count);
        mexPrintf("Actual number of points: %d.\n",Num_Points);
        mexErrMsgTxt("Number of points printed out does not match the actual number of points!\n");
        }*/
}

inline double TR::compare(const d_node *p1, const d_node* p2) 
{
  return p1 -> dist - p2 -> dist;
}

void TR::halfsort (v_array<d_node > cover_set)
{
  if (cover_set.index <= 1)
    return;
  register d_node *base_ptr =  cover_set.elements;

  d_node *hi = &base_ptr[cover_set.index - 1];
  d_node *right_ptr = hi;
  d_node *left_ptr;
  
  while (right_ptr > base_ptr)
    {
      d_node *mid = base_ptr + ((hi - base_ptr) >> 1);
      
      if (compare ( mid,  base_ptr) < 0.)
	SWAP (mid, base_ptr);
      if (compare ( hi,  mid) < 0.)
	SWAP (mid, hi);
      else
	goto jump_over;
      if (compare ( mid,  base_ptr) < 0.)
	SWAP (mid, base_ptr);
    jump_over:;
      
      left_ptr  = base_ptr + 1;
      right_ptr = hi - 1;
      
      do
	{
	  while (compare (left_ptr, mid) < 0.)
	    left_ptr ++;
	  
	  while (compare (mid, right_ptr) < 0.)
	    right_ptr --;
	  
	  if (left_ptr < right_ptr)
	    {
	      SWAP (left_ptr, right_ptr);
	      if (mid == left_ptr)
		mid = right_ptr;
	      else if (mid == right_ptr)
		mid = left_ptr;
	      left_ptr ++;
	      right_ptr --;
	    }
	  else if (left_ptr == right_ptr)
	    {
	      left_ptr ++;
	      right_ptr --;
	      break;
	    }
	}
      while (left_ptr <= right_ptr);
      
      hi = right_ptr;
    }
}

v_array<v_array<d_node> > TR::get_cover_sets(v_array<v_array<v_array<d_node> > > &spare_cover_sets)
{
  v_array<v_array<d_node> > ret = pop(spare_cover_sets);
  while (ret.index < 151)
    {
      v_array<d_node> temp;
      push(ret, temp);
    }
  return ret;
}

inline bool TR::shell(double parent_query_dist, double child_parent_dist, double upper_bound)
{
  return parent_query_dist - child_parent_dist <= upper_bound;
  //    && child_parent_dist - parent_query_dist <= upper_bound;
}

void TR::update(double *k_upper_bound, double upper_bound)
{
  double *end = k_upper_bound;// + internal_k-1;
  double *begin = k_upper_bound;
  for (;end != begin; begin++)
    {
      if (upper_bound < *(begin+1))
	*begin = *(begin+1);
      else {
	*begin = upper_bound;
	break;
      }
    }
  if (end == begin)
    *begin = upper_bound;
}
double *TR::alloc_upper()
{
  return (double *)malloc(sizeof(double));
}
void TR::setter(double* begin, double max)
{
    *begin = max;
}

//void (*update)TR::(double *foo, double bar) = update_k;
//void (*setter)TR::(double *foo, double bar) = set_k;
//double* (*alloc_upper)TR::() = alloc_epsilon_approx;

inline void TR::copy_zero_set(node* query_chi, double* new_upper_bound, 
			    v_array<d_node> &zero_set, v_array<d_node> &new_zero_set)
{
  new_zero_set.index = 0;
  d_node *end = zero_set.elements + zero_set.index;
  for (d_node *ele = zero_set.elements; ele != end ; ele++)
    {
      double upper_dist = *new_upper_bound + 2. * query_chi->max_dist;
      if (shell(ele->dist, query_chi->parent_dist, upper_dist))
	{
	  double d = Point_Obj->distance(query_chi->p, ele->n->p, upper_dist,&counter);
	  
	  if (d <= upper_dist)
	    {
	      if (d < *new_upper_bound) 
		update(new_upper_bound, d);
	      d_node temp = {d, ele->n};
	      push(new_zero_set,temp);
	    }
	}
    }
}

inline void TR::copy_cover_sets(node* query_chi, double* new_upper_bound,
			      v_array<v_array<d_node> > &cover_sets,
			      v_array<v_array<d_node> > &new_cover_sets,
			      int current_scale, int max_scale)
{
  for (; current_scale <= max_scale; current_scale++)
    {
      d_node* ele = cover_sets[current_scale].elements;
      d_node* end = cover_sets[current_scale].elements + cover_sets[current_scale].index;
      for (; ele != end; ele++)
	{ 
	  double upper_dist = *new_upper_bound + 2. *query_chi->max_dist + ele->n->max_dist;
	  if (shell(ele->dist, query_chi->parent_dist, upper_dist))
	    {
	      double d = Point_Obj->distance(query_chi->p, ele->n->p, upper_dist, &counter);
	      
	      if (d <= upper_dist)
		{
		  if (d < *new_upper_bound)
		    update(new_upper_bound,d);
		  d_node temp = {d, ele->n};
		  push(new_cover_sets[current_scale],temp);
		}
	    }
	}
    }  
}

void TR::print_query(const node *top_node)
{
  mexPrintf ("query = \n");
  point_str *p = top_node->p;
  Point_Obj->print(p);
  if ( top_node->num_children > 0 ) {
    mexPrintf("scale = %i\n",top_node->scale);
    mexPrintf("max_dist = %f\n",top_node->max_dist);
    mexPrintf("num children = %i\n",top_node->num_children);
  }
}

void TR::print_cover_sets(v_array<v_array<d_node> > &cover_sets,
		      v_array<d_node> &zero_set,
		      int current_scale, int max_scale)
{
  mexPrintf("cover set = \n");
  for (; current_scale <= max_scale; current_scale++)
    {
      d_node* ele = cover_sets[current_scale].elements;
      d_node* end = cover_sets[current_scale].elements + cover_sets[current_scale].index;
      mexPrintf ("%i\n", current_scale);
      for (; ele != end; ele++)
	{
	  node *n = (node *)ele->n;
	  Point_Obj->print(n->p);
	}
    }
  d_node *end = zero_set.elements + zero_set.index;
  mexPrintf ("infinity\n");
  for (d_node *ele = zero_set.elements; ele != end ; ele++)
    {
      node *n = (node *)ele->n;
      Point_Obj->print(n->p);
    }
}


/*
  An optimization to consider:
  Make all distance evaluations occur in descend.

  Instead of passing a cover_set, pass a stack of cover sets.  The
  last element holds d_nodes with your distance.  The next lower
  element holds a d_node with the distance to your query parent,
  next = query grand parent, etc..

  Compute distances in the presence of the tighter upper bound.
 */

inline void TR::descend(const node* query, double* upper_bound, 
		      int current_scale,
		      int &max_scale, v_array<v_array<d_node> > &cover_sets, 
		      v_array<d_node> &zero_set)
{
  d_node *end = cover_sets[current_scale].elements + cover_sets[current_scale].index;
  for (d_node *parent = cover_sets[current_scale].elements; parent != end; parent++)
    {
      const node *par = parent->n;
      double upper_dist = *upper_bound + query->max_dist + query->max_dist;
      if (parent->dist <= upper_dist + par->max_dist && par->children)
	{
	  node *chi = par->children;
	  if (parent->dist <= upper_dist + chi->max_dist)
	    {
	      if (chi->num_children > 0)
		{
		  if (max_scale < chi->scale)
		    max_scale = chi->scale;
		  d_node temp = {parent->dist, chi};
		  push(cover_sets[chi->scale], temp);
		}
	      else if (parent->dist <= upper_dist)
		{
		  d_node temp = {parent->dist, chi};
		  push(zero_set, temp);
		}
	    }
	  node *child_end = par->children + par->num_children;
	  for (chi++; chi != child_end; chi++)
	    {
	      double upper_chi = *upper_bound + chi->max_dist + query->max_dist + query->max_dist;
	      if (shell(parent->dist, chi->parent_dist, upper_chi))
		{
		  double d = Point_Obj->distance(query->p, chi->p, upper_chi, &counter);
		  if (d <= upper_chi) 
		    {
		      if (d < *upper_bound)
			update(upper_bound, d);
		      if (chi->num_children > 0)
			{
			  if (max_scale < chi->scale)
			    max_scale = chi->scale;
			  d_node temp = {d, chi};
			  push(cover_sets[chi->scale],temp);
			}
		      else 
			if (d <= upper_chi - chi->max_dist)
			  {
			    d_node temp = {d, chi};
			    push(zero_set, temp);
			  }
		    }
		}
	    }
	}
    }
}


inline void TR::descend_approx(const node* query, double* upper_bound, 
		      int current_scale,
		      int &max_scale, v_array<v_array<d_node> > &cover_sets, 
		      v_array<d_node> &zero_set)
{
  d_node *end = cover_sets[current_scale].elements + cover_sets[current_scale].index;
  for (d_node *parent = cover_sets[current_scale].elements; parent != end; parent++)
    {
      const node *par = parent->n;
      if(internal_epsilon_approx>0 && query->scale == 150)
      {    
          if (parent->dist >= (par->max_dist + query->max_dist + query->max_dist)*(1+1/internal_epsilon_approx))
          {
	          d_node temp = {parent->dist, par};
              push(zero_set, temp);
              return;
          }
      }
      if(stop_scale>0)
      {    
          if (par->scale >= stop_scale)
          {
	          d_node temp = {parent->dist, par};
              push(zero_set, temp);
              return;
          }
      }
      double upper_dist = *upper_bound + query->max_dist + query->max_dist;
      if (parent->dist <= upper_dist + par->max_dist && par->children)
	{
	  node *chi = par->children;
	  if (parent->dist <= upper_dist + chi->max_dist)
	    {
	      if (chi->num_children > 0)
		{
		  if (max_scale < chi->scale)
		    max_scale = chi->scale;
		  d_node temp = {parent->dist, chi};
		  push(cover_sets[chi->scale], temp);
		}
	      else if (parent->dist <= upper_dist)
		{
		  d_node temp = {parent->dist, chi};
		  push(zero_set, temp);
		}
	    }
	  node *child_end = par->children + par->num_children;
	  for (chi++; chi != child_end; chi++)
	    {
	      double upper_chi = *upper_bound + chi->max_dist + query->max_dist + query->max_dist;
	      if (shell(parent->dist, chi->parent_dist, upper_chi))
		{
		  double d = Point_Obj->distance(query->p, chi->p, upper_chi,&counter);
		  if (d <= upper_chi) 
		    {
		      if (d < *upper_bound)
			update(upper_bound, d);
		      if (chi->num_children > 0)
			{
			  if (max_scale < chi->scale)
			    max_scale = chi->scale;
			  d_node temp = {d, chi};
			  push(cover_sets[chi->scale],temp);
			}
		      else 
			if (d <= upper_chi - chi->max_dist)
			  {
			    d_node temp = {d, chi};
			    push(zero_set, temp);
			  }
		    }
		}
	    }
	}
    }
}

void TR::brute_nearest(const node* query,v_array<d_node> zero_set,
		   double* upper_bound,
		   v_array<v_array<point> > &results, v_array<double> &dist_mins,
		   v_array<v_array<d_node> > &spare_zero_sets)
{
  if (query->num_children > 0)
    {
      v_array<d_node> new_zero_set = pop(spare_zero_sets);
      node* query_chi = query->children; 
      brute_nearest(query_chi, zero_set, upper_bound, results, dist_mins, spare_zero_sets);
      double* new_upper_bound = alloc_upper();
      
      node *child_end = query->children + query->num_children;
      for (query_chi++;query_chi != child_end; query_chi++)
	{
	  setter(new_upper_bound,*upper_bound + query_chi->parent_dist);
	  copy_zero_set(query_chi, new_upper_bound, zero_set, new_zero_set);
	  brute_nearest(query_chi, new_zero_set, new_upper_bound, results, dist_mins, spare_zero_sets);
	}
      free (new_upper_bound);
      new_zero_set.index = 0;
      push(spare_zero_sets, new_zero_set);
    }
  else 
    {
      v_array<point> temp;
      push(temp, query->p);
      d_node *end = zero_set.elements + zero_set.index;
      for (d_node *ele = zero_set.elements; ele != end ; ele++)
      {
	      if (ele->dist <= *upper_bound) 
	      {
              push(temp, ele->n->p);
              push(dist_mins, upper_bound[0]);//ele->dist);
              break;//added, so only one nearest neighor will be found...
          }
      }
      push(results,temp);
      //push(dmins,upper_bound[0]);
    }
}

void TR::internal_batch_nearest_neighbor(const node *query, 
				     v_array<v_array<d_node> > &cover_sets,
				     v_array<d_node> &zero_set,
				     int current_scale,
				     int max_scale,
				     double* upper_bound,
				     v_array<v_array<point> > &results,
                     v_array<double> &dist_mins,
				     v_array<v_array<v_array<d_node> > > &spare_cover_sets,
				     v_array<v_array<d_node> > &spare_zero_sets)
{
  if (current_scale > max_scale) // All remaining points are in the zero set. 
    brute_nearest(query, zero_set, upper_bound, results, dist_mins, spare_zero_sets);
  else
    if (query->scale <= current_scale && query->scale != 150) 
      // Our query has too much scale.  Reduce.
      { 
	node *query_chi = query->children;
	v_array<d_node> new_zero_set = pop(spare_zero_sets);
	v_array<v_array<d_node> > new_cover_sets = get_cover_sets(spare_cover_sets);
	double* new_upper_bound = alloc_upper();

	node *child_end = query->children + query->num_children;
	for (query_chi++; query_chi != child_end; query_chi++)
	  {
	    setter(new_upper_bound,*upper_bound + query_chi->parent_dist);
	    copy_zero_set(query_chi, new_upper_bound, zero_set, new_zero_set);
	    copy_cover_sets(query_chi, new_upper_bound, cover_sets, new_cover_sets,
			      current_scale, max_scale);
	    internal_batch_nearest_neighbor(query_chi, new_cover_sets, new_zero_set,
					    current_scale, max_scale, new_upper_bound, 
					    results, dist_mins, spare_cover_sets, spare_zero_sets);
	  }
	free (new_upper_bound);
	new_zero_set.index = 0;
	push(spare_zero_sets, new_zero_set);
	push(spare_cover_sets, new_cover_sets);
	internal_batch_nearest_neighbor(query->children, cover_sets, zero_set, 
					current_scale, max_scale, upper_bound, results, dist_mins,
					spare_cover_sets, spare_zero_sets);
      }
   /* else if (query->scale > current_scale && query->scale != 150)
    {
 	halfsort(cover_sets[current_scale]);
	descend(query, upper_bound, current_scale, max_scale,cover_sets, zero_set);
	cover_sets[current_scale++].index = 0;
	internal_batch_nearest_neighbor(query, cover_sets, zero_set, 
					current_scale, max_scale, upper_bound, results, dist_mins,
					spare_cover_sets, spare_zero_sets);
    }*/
    else // reduce cover set scale
      {
	halfsort(cover_sets[current_scale]);
	descend_approx(query, upper_bound, current_scale, max_scale,cover_sets, zero_set);
	cover_sets[current_scale++].index = 0;
	internal_batch_nearest_neighbor(query, cover_sets, zero_set, 
					current_scale, max_scale, upper_bound, results, dist_mins,
					spare_cover_sets, spare_zero_sets);
      }
}

void TR::batch_nearest_neighbor(const node &top_node, const node &query, 
			    v_array<v_array<point> > &results, v_array<double> &dist_mins)
{
  v_array<v_array<v_array<d_node> > > spare_cover_sets;
  v_array<v_array<d_node> > spare_zero_sets;

  v_array<v_array<d_node> > cover_sets = get_cover_sets(spare_cover_sets);
  v_array<d_node> zero_set = pop(spare_zero_sets);
  
  double* upper_bound = alloc_upper();
  setter(upper_bound,MAXFLOAT);
  
  double top_dist = Point_Obj->distance(query.p, top_node.p, MAXFLOAT,&counter);
  update(upper_bound, top_dist);
  d_node temp = {top_dist, &top_node};
  push(cover_sets[0], temp);
  
  internal_batch_nearest_neighbor(&query,cover_sets,zero_set,0,0,upper_bound,results, dist_mins,
				  spare_cover_sets,spare_zero_sets);
  
  free(upper_bound);
  push(spare_cover_sets, cover_sets);
  
  for (int i = 0; i < spare_cover_sets.index; i++)
    {
      v_array<v_array<d_node> > cover_sets = spare_cover_sets[i];
      for (int j = 0; j < cover_sets.index; j++)
	free (cover_sets[j].elements);
      free(cover_sets.elements);
    }
  free(spare_cover_sets.elements);
  
  push(spare_zero_sets, zero_set);

  for (int i = 0; i < spare_zero_sets.index; i++)
    free(spare_zero_sets[i].elements);
  free(spare_zero_sets.elements);
}

void TR::batch_nearest_neighbor_current(const node &top_node, const node &query, 
			    v_array<v_array<point> > &results, v_array<double> &dist_mins, double input_uppbound)
{
  v_array<v_array<v_array<d_node> > > spare_cover_sets;
  v_array<v_array<d_node> > spare_zero_sets;

  v_array<v_array<d_node> > cover_sets = get_cover_sets(spare_cover_sets);
  v_array<d_node> zero_set = pop(spare_zero_sets);
  
  double* upper_bound = alloc_upper();
  setter(upper_bound,input_uppbound);
  
  double top_dist = Point_Obj->distance(query.p, top_node.p, MAXFLOAT,&counter);
  update(upper_bound, top_dist);
  d_node temp = {top_dist, &top_node};
  push(cover_sets[0], temp);
  
  internal_batch_nearest_neighbor(&query,cover_sets,zero_set,0,0,upper_bound,results, dist_mins,
				  spare_cover_sets,spare_zero_sets);

  //dmin = upper_bound[0];
  
  free(upper_bound);
  push(spare_cover_sets, cover_sets);
  
  for (int i = 0; i < spare_cover_sets.index; i++)
    {
      v_array<v_array<d_node> > cover_sets = spare_cover_sets[i];
      for (int j = 0; j < cover_sets.index; j++)
	free (cover_sets[j].elements);
      free(cover_sets.elements);
    }
  free(spare_cover_sets.elements);
  
  push(spare_zero_sets, zero_set);

  for (int i = 0; i < spare_zero_sets.index; i++)
    free(spare_zero_sets[i].elements);
  free(spare_zero_sets.elements);
}

/*void TR::epsilon_approx_nearest_neighbor(const node &top_node, const node &query, 
			      v_array<v_array<point> > &results, double epsilon_approx)
{
  internal_epsilon_approx = epsilon_approx;
  update = update_epsilon_approx;
  setter = set_epsilon_approx;
  alloc_upper = alloc_epsilon_approx;

  batch_nearest_neighbor(top_node, query,results);
}*/
/***************************************************************************************/



#undef  TR

/***/
