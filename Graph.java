import java.util.LinkedList;
import java.util.ArrayList;
import java.util.Scanner;

import java.util.Queue;

import java.util.Arrays;

import java.io.File;
import java.io.IOException;

public class Graph
{
	// We want to have Topo Sort for a graph
	// How this works, is such.
	// Each node will have a "rank" of dependencies, indicating how many Nodes are required to reach it.
	// Every time we encounter a node that has 0 dependencies, we will lower the dependency value of each Node it is adjacent to
	public void topologicalSort()
	{
		boolean [] used = new boolean[n];

		int [] path = new int[n];
		int [] dependencies = new int[n];

		// Fill dependency list will amount of edges incident to it
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (matrix[j][i] == true)
				{
					dependencies[i]++;
				}
			}
		}

		topologicalSort(used, dependencies, path, 0);
	}

	private boolean topologicalSort(boolean [] used, int [] dependencies, int [] path, int index)
	{
		if (index == n)
		{
			printPath(path);
			return true;
		}

		// For every vertex in the graph, let's check its dependecy values
		for (int i = 0; i < n; i++)
		{
			// if we have no dependencies and current node is not used
			if (dependencies[i] == 0 && !used[i])
			{
				// Lower the dependecy value for each Node adjacent to i'th node
				for (int j = 0; j < n; j++)
				{
					if (matrix[i][j] == true)
					{
						dependencies[j]--;
					}
				}
			}

			// Set this Node as used, update path to print, and run backtrack state change
			used[i] = true;

			path[index] = i;

			if (topologicalSort(used, dependencies, path, index + 1))
			{
				return true;
			}

			// Revert changes
			used[i] = false;

			for (int j = 0; j < n; j++)
			{
				if (matrix[i][j] == true)
				{
						dependencies[j]++;
				}
			}
		}

		return false;
	}

	// This works in a very similar way to the other topo sort implementation. Except now we are using a List.
	public void topologicalSortWList()
	{
		boolean [] used = new boolean[n];
		
		int [] path = new int[n];
		int [] dependencies = new int[n];

		for (ArrayList<Integer> adj : list)
		{
			for (Integer node : adj)
			{
				dependencies[node]++;
			}
		}

		topologicalSortWList(used, dependencies, path, 0);
	}

	// Same as other topo sort implementation except utilizing arraylist
	private boolean topologicalSortWList(boolean [] used, int [] dependencies, int [] path, int index)
	{
		if (index == n)
		{
			printPath(path);
			return true;
		}

		for (int i = 0; i < n; i++)
		{
			if (dependencies[i] == 0 && !used[i])
			{
				for (Integer node : list.get(i))
				{
					dependencies[node]--;
				}

				used[i] = true;
				path[index] = i;

				if (topologicalSortWList(used, dependencies, path, index+1))
				{
					return true;
				}

				used[i] = false;

				for (Integer node : list.get(i))
				{
					dependencies[node]++;
				}
			}
		}

		return false;
	}
			
	
	// We want to have BFS for graph recursively
        public void bfs(int node)
        {
            boolean [] used = new boolean[n];

            int [] path = new int[n];

            used[node] = true;
            path[0] = node;

            bfs(used, path, node, 1);
        }

        // In here, we want to mimick the queue operation by having a Set of Nodes to visit recursive call.
        // Since we are already doing recursion, might as well print all possible combinations as well :)
        //
        // To print all possible BFS's, we have to permutate each "Children" ArrayList so that we get different orderings.
        // How do we permute elements?
        //
        // 
        
        private void bfs(boolean [] used, int [] path, int node, int index)
        {
            // If we tracked all elements
            //
            if (index == n)
            {
                printPath(path);
                return;
            }

            // For each children of current state, we want to create a List of permuted elements for each.
            // So for example, we are at Node 0 (root). Root has 3 children. Then we can simply create all permutations of 
            // those children. But once we BFS'd down to where we are finding the children of all those BF's,
            // we need to find the permutations for each Node's children seperately and join them according to the
            // original order of it's parent.
            
            // This arraylist holds all the possible permutations for each Node's children

            // We first want to find all the children of current Node
            ArrayList<Integer> list = new ArrayList<>();

            // For all vertices, we want to find which are adjacent to current Node
            for (int i = 0; i < n; i++)
            {
                // We only want to add Nodes that have not been used
                // The i'th term indicates a Node that we are adjacent to
                if (matrix[node][i] && !used[i])
                {
                    list.add(i);
                }
            }
            
            // Create arraylist that will hold all the permutations of current Node's children
            ArrayList<ArrayList<Integer>> perms = new ArrayList<>();

            // Permute the list, which is stored in perms
            permuteList(perms, list, list.size());

            // We are iterating through each permutation
            for (ArrayList<Integer> l : perms)
            {
                // Through each specfic node we can visit in that permutation
                for (Integer nextNode : l)
                {
                    // We are visiting that node
                    used[nextNode] = true;
                    path[index] = nextNode;
                    bfs(used, path, nextNode, index + 1);

                    used[nextNode] = false;
                }
            }
        }
	
	// We want to have BFS iteratively
        public void bfsIteratively(int start)
        {
            boolean [] visited = new boolean[n];

            Queue<Integer> q = new LinkedList<>();

            // we want to add to the queue, and for every pop, we want to get all the adjacent nodes to it

            visited[start] = true;
            q.add(start);

            // we want to first visit all the Nodes for current Node
            // we do that by popping from queue
            while (!q.isEmpty())
            {
                int node = q.remove();
                printNode(node);

                // we want to add all of the Current Node's children to visit first
                for (int i = 0; i < n; i++)
                {
                    // if we havent visited this node yet
                    if (matrix[node][i] && !visited[i])
                    {
                        visited[i] = true;
                        q.add(i);
                    }
                }
            }
        }

        // Add all permutations of a Set of numbers to a List, includes original List
        private void permuteList(ArrayList<ArrayList<Integer>> permuted, 
                                 ArrayList<Integer> list, int size)
        {
            // Heap's base case
            if (size == 1)
            {
                permuted.add(new ArrayList<Integer>(list));
                return;
            }

            // Along each step of the way, we are permuting from left-to-right
            for (int i = 0; i < size; i++)
            {
                // we permute this state, and the state where we permute the element next to this one
                // so we permute [ab]c
                // and a[bc]
                //
                // I think... probably wrong tho
                permuteList(permuted, list, size - 1);
                
                // if even, we switch first and "last" element, indicated by size. Realize last is
                // not the last element of the list.
                if (size%2==0)
                {
                    int temp = list.get(i);
                    list.set(i, list.get(size-1));
                    list.set(size-1, temp);
                }

                // if odd, we switch first and "last" element
                else
                {
                    int temp = list.get(0);
                    list.set(0, list.get(size-1));
                    list.set(size-1, temp);
                }
            }
        }
        
        /*
        private void swap(Integer a, Integer b)
        {
            int temp = a;
            a = b;
            b = temp;
        }
        */
    
        public void printNode(int node)
        {
            System.out.print(value[node]);
        }

	
	// We want to have DFS for graph recursively
	
	// We want to have DFS iteratively
	
	// Check if graph is bipartite
	
	// Find graph is k-colorable
	
	// Count number of connected componenets
	
	// Dijkstras implementation
	
	// Bellman-ford
	
	// Kruskals algo
	
	// Prims algo
	
	// Find Hamiltonian Path / Cycle
	
	// Find Eulerian Path / Cycle ?

	// Do the same for all, with Lists! Good luck
	
	// MORE NOTES:
	// Try to follow along more of the stuff on Webcourses. Here are some things to try.
	// For every algorithm implementation. Try to get a single solution, solution; and a multi-solution, solution.
	// Modify my implementation of Dijkstra's algorithm to print its output in alphabetic order by vertex.
	// Modify my implementation of Dijkstra's algorithm to print not just the shortest path lengths, but also the actual paths taken.
	//
	// Code up an O(|V|2) version of Dijkstra's algorithm that does not use a minheap. This might be trickier than you'd think; you have to 
	// use an array where removing the smallest element happens in O(1) time (after finding that smallest element in O(n) time). 
	// So, you can't just use an ArrayList and call the remove() method to pull out the smallest element from an arbitrary position, 
	// since removing the first element of an ArrayList is an O(n) operation
	//
	// Use Bellman-Ford to track negative cycles.
	// Maybe? Code up single-room scheduling problem, multiple-room scheduling problem, continuous knapsack problem.
	//
	// Backtracking side:
	// Fox-Goose-Problem solve. (Already solved it)
	// Missionaries and Cannibals Problem.
	// Iterative solution to NQueens problem. (Try to also get it print only one solution, many solutions, and only two solutions without class member variables) - Maybe with recursion and without ;)
	// Follow along 7/8 in #Backtracking
	//
	
	// Constructor that loads a map from file
	
	private void printPath(int [] path)
	{
		for (int i = 0; i < path.length - 1; i++)
		{
			System.out.print(value[i] + " -> ");
		}
		
		if (path.length > 1)
		{
			System.out.println(value[path.length-1]);
		}
	}

	int n;
	boolean [][] matrix;

	String [] value;
	ArrayList<ArrayList<Integer>> list;

	Graph(String file) throws IOException
	{
		Scanner scan = new Scanner(new File(file));
		
		n = scan.nextInt();

		// Create list and matrix array for graph presentation
		list = new ArrayList<>();
		matrix = new boolean[n][n];
		value = new String[n];
		
		for (int i = 0; i < n; i++)
		{
			// This is a sublist for one of the Vertices, indicating all the nodes adjacent to it
			ArrayList<Integer> subList = new ArrayList<>();
			
			value[i] = scan.next();

			// Update matrix as needed
			for (int j = 0; j < n; j++)
			{
				int val = scan.nextInt();
				matrix[i][j] = (val == 1);

				if (val == 1)
				{
					subList.add(j);
				}
			}

			list.add(subList);
		}
	}

	public static void main(String [] args) throws IOException
	{
		Graph g = new Graph(args[0]);
		g.topologicalSort();
		g.topologicalSortWList();
                g.bfsIteratively(0);
                g.bfs(0);
        }
}
