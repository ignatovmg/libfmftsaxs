#include <stdio.h>
#include <malloc.h>
#include "ptree.h"

/*void printProfileImage(FILE* stream, profileImage* im)
{
	fprintf(stream, "\nDistance: %.1f\n1st protein rotated: (0.0, %.1f, %.1f)\n2nd protein rotated: (%.1f, %.1f, %.1f)\nScore = %.5f\n",
		im->u, 
		im->theta, 
		im->fi,
		im->alpha, 
		im->beta, 
		im->gamma,
		im->score);
}*/

int compareImgByScore(profileImage* a, profileImage* b) 
{
	return (int)(a->score*10000.0 - b->score*10000.0);
}

tree* treeCreate(void)
{
	tree* new_tree = (tree*)malloc(sizeof(tree));
	if (new_tree == NULL) return NULL;
	new_tree->root = NULL;
	new_tree->count = 0;
	return new_tree;
}

int insert(tree* search_tree, profileImage* item)
{
	node* search_node;
	node** nw;

	nw = &search_tree->root;
	search_node = search_tree->root;

	while(1)
	{
		if(search_node == NULL)
		{
			search_node = *nw = (node*)malloc(sizeof(search_node));
			if(search_node != NULL)
			{
				search_node->img = item;
				search_node->left = search_node->right = NULL;
				search_tree->count++;
				return 1;
			}
			else 
				return 0;
		}
		else 
			if(compareImgByScore(item, search_node->img) >= 0)
			{
				nw = &search_node->right;
				search_node = search_node->right;
			}
			else
			{
				nw = &search_node->left;
				search_node = search_node->left;
			}
	}
}

void deleteLast(tree* search_tree)
{
	node* search_node = search_tree->root;
	node** link = &search_tree->root;

	if(search_node == NULL) 
		return;
	if (search_node->right == NULL)
	{
		node* p = *link;
		*link = (*link)->left;
		free(p->img);
		free(p);
	}
	else
	{
		node* p;
		while(search_node->right != NULL)
		{
			p = search_node;
			search_node = search_node->right;
		}
		p->right = search_node->left;
		free(search_node->img);
		free(search_node);
	}
	search_tree->count--;
}

static void walk(const node* search_node)
{
	if(search_node == NULL) 
		return;
	walk(search_node->left);
	printProfileImage(stdout, search_node->img);
	walk(search_node->right);
}	

void walkAll(const tree* my_tree)
{
	walk(my_tree->root);
}

static void deallocateNode(node* n)
{
	if (n->img != NULL)
		free(n->img);
	if (n->left != NULL)
		deallocateNode(n->left);
	if (n->right != NULL)
		deallocateNode(n->right);
	free(n);
	return;
}

void deallocateTree(tree* search_tree)
{
	deallocateNode(search_tree->root);
	free(search_tree);
}

static void untree(node* search_node, profileImage** ar, int* index)
{
	if(search_node == NULL) 
		return;
	untree(search_node->left, ar, index);
	ar[*index] = search_node->img;
	(*index)++;
	untree(search_node->right, ar, index);
}	

profileImage** untreeAll(tree* my_tree, profileImage** ar)
{
	int index = 0;
	untree(my_tree->root, ar, &index);
	return ar;
}

/*void main(void)
{
	tree* tr = treeCreate();
	//profileImage* ims[10];
	int j;
	for (j = 0; j < 3; j++)
	{
		profileImage* i = malloc(sizeof(profileImage)); i->score = j;
		printf("j = %i\n", j);
		insert(tr, i);
	}
	//walkAll(tr);
	printf("jlkjhlkj\n");
	deleteLast(tr);
	printf("jlkjhlkj\n");
	deleteLast(tr);
	printf("jlkjhlkj\n");
	profileImage* i = malloc(sizeof(profileImage)); i->score = 4; i->alpha = 20;
		insert(tr, i);
	i = malloc(sizeof(profileImage)); i->score = 2; i->alpha = 20; insert(tr, i);
	deleteLast(tr);
	i = malloc(sizeof(profileImage)); i->score = -1; i->alpha = 20; insert(tr, i);
	i = malloc(sizeof(profileImage)); i->score = -2; i->alpha = 20; insert(tr, i);
	deleteLast(tr);
	//deleteLast(tr);
	//deleteLast(tr);
	printf("jlkjhlkj\n");
	printf("===============\n");
	walkAll(tr);
	printf("===============\n");
	profileImage** ar = calloc(tr->count, sizeof(profileImage*));
	untreeAll(tr, ar);
	printf("count = %i\n", tr->count);
	for (j = 0; j < tr->count; j++)
		printProfileImage(stdout, ar[j]);
}*/

