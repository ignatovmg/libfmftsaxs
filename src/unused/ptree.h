#ifndef _PTREE_H
#define _PTREE_H

/**
* API for BST with compressed MasshaProfiles (profileImage)
*/

typedef struct profileImage_
{
	double score;
	double u;
	int alpha;
	int beta;
	int gamma;
	int theta;
	int fi;
}profileImage;

void printProfileImage(FILE* stream, profileImage* im);
int compareImgByScore(profileImage* a, profileImage* b);

typedef struct node_
{
	profileImage* img;
	
	struct node_* left;
	struct node_* right;
}node;

typedef struct tree_
{
	node* root;	
	int count;			
}tree;

tree* treeCreate(void);
int insert(tree* search_tree, profileImage* item);
void deleteLast(tree* search_tree);
//static void walk(const node* search_node);
void walkAll(const tree* my_tree);
//static void deallocateNode(node* n);
void deallocateTree(tree* search_tree);
//static void untree(node* search_node, profileImage** ar, int* index);
profileImage** untreeAll(tree* my_tree, profileImage** ar);
#endif
