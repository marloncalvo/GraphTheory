import java.util.LinkedList;
import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.Queue;
import java.util.Stack;
import java.util.Arrays;

import java.util.Scanner;

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
        System.out.println("All Topo Sort Traversals...");
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

    public void topologicalSortIteratively(int startNode)
    {
        int [] dependencies = new int[n];
        int [] path = new int[n];

        Stack<Integer> s = new Stack<>();

        for (ArrayList<Integer> adj : list)
        {
            for (Integer adjNode : adj)
            {
                dependencies[adjNode]++;
            }
        }
        System.out.println("Topo Sort Iteratively...");

        s.push(startNode);

        int i = 0;
        while (!s.isEmpty())
        {
            int node = s.pop();
            path[i++] = node;
            // want to deincrement priority for all nodes adjecnt to it
            for (Integer adjNode : list.get(node))
            {
                dependencies[adjNode]--;
                if (dependencies[adjNode] == 0)
                {
                    s.push(adjNode);
                }
            }
        }

        System.out.println("\tFound \"path\"...");
        System.out.print("\t");
        printPath(path);
        System.out.println();
    }

	private void topologicalSort(boolean [] used, int [] dependencies, int [] path, int index)
	{
		if (index == n)
		{
            System.out.println("\tFound \"path\"...");
            System.out.print("\t");
			printPath(path);
            System.out.println();
			return;
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

    			// Set this Node as used, update path to print, and run backtrack state change
    			used[i] = true;

    			path[index] = i;

    			topologicalSort(used, dependencies, path, index + 1);

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
		}
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
        System.out.println("Topo Sort Traversal...");

		topologicalSortWList(used, dependencies, path, 0);
	}

	// Same as other topo sort implementation except utilizing arraylist
	private boolean topologicalSortWList(boolean [] used, int [] dependencies, int [] path, int index)
	{
		if (index == n)
		{
			System.out.println("\tFound path...");
            System.out.print("\t");
            printPath(path);
            System.out.println();
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
	
	// We want to have BFS iteratively
        public void bfsIteratively(int start)
        {
            boolean [] visited = new boolean[n];

            Queue<Integer> q = new LinkedList<>();

            int [] path = new int[n];

            // we want to add to the queue, and for every pop, we want to get all the adjacent nodes to it

            System.out.println("BFS Iterative Traversal...");

            visited[start] = true;
            q.add(start);

            // we want to first visit all the Nodes for current Node
            // we do that by popping from queue
            int k = 0;
            while (!q.isEmpty())
            {
                int node = q.remove();
                path[k++] = node;

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

            System.out.println("\tFound path...");
            System.out.print("\t");
            printPath(path);
            System.out.println();
        }

        // We want to have BFS for graph recursively
        public void bfs(int node)
        {
            boolean [] used = new boolean[n];
            
            int [] path = new int[n];
            Arrays.fill(path, -1);

            used[node] = true;
            path[0] = node;

            ArrayList<Integer> nodes = new ArrayList<>();
            nodes.add(0);

            System.out.println("All BFS Traversals...");
            bfs(used, path, nodes, 1);
        }

        // In here, we want to mimick the queue operation by having a Set of Nodes to visit recursive call.
        // Since we are already doing recursion, might as well print all possible combinations as well :)
        //
        // To print all possible BFS's, we have to permutate each "Children" ArrayList so that we get different orderings.
        // How do we permute elements?
        //
        // 
        
        // This method prints all BFS's possible with recursion.
        // We must make sure that we are printing all the permutations possible for each Node's children.
        // What we are doing is from the given list of Nodes, we will add those to the path first.
        // Then, we want to visit a Node, create a permutation of it, and then visit all other Nodes and visit all possible permutations of each.
        // Visit those paths. And do the same for each permutation of that 1st node.
        private void bfs(boolean [] used, int [] path, ArrayList<Integer> nodes, int index)
        {
            // If we tracked all elements
            if (index >= n)
            {
                System.out.println("\tFound path...");
                System.out.print("\t");
                printPath(path);
                System.out.println("");
                return;
            }

            // Store all the children for each Node seperately
            ArrayList<ArrayList<Integer>> children = new ArrayList<>();


            //System.out.println("INFO FOR CURRENT STATE OF BFS\t");
            //printPathAsSet(nodes);
            //System.out.println();
            //System.out.println("\tCHILDREN FOR NODE");

            // For all vertices, we want to find which are adjacent to current Node
            for (Integer node : nodes)
            {
                //System.out.print("\t");
                //printNode(node);
                ArrayList<Integer> children_t = new ArrayList<>();
                for (int i = 0; i < n; i++)
                {
                    // We only want to add Nodes that have not been used
                    // The i'th term indicates a Node that we are adjacent to
                    if (matrix[node][i] && !used[i])
                    {
                        children_t.add(i);
                    }
                }

                //System.out.print("\t");
                //printPathAsSet(children_t);
                children.add(children_t);
            }

            //System.out.println();

            // We are iterating through each permutation
            for (ArrayList<Integer> perms : permuteNodes(children))
            {
                /*
                System.out.println("\tPermutation...");
                System.out.print("\t");
                printPathAsSet(perms);
                System.out.println();
                */

                // Through each specfic node we can visit in that permutation
                int i = 0;

                //System.out.println(index);
                for (Integer nextNode : perms)
                {
                    // We are visiting that node
                    used[nextNode] = true;
                    path[index+i] = nextNode;
                    i++;
                }
/*
                System.out.println("PRINTING CURRENT PATH");
                System.out.print("\t");
                printPath(path);
                System.out.println();
*/
                bfs(used, path, perms, index + i);

                for (Integer nextNode : perms)
                {
                    used[nextNode] = false;
                }
            }
        }

        // We want to store all the possible permutation for each Node
        private ArrayList<ArrayList<Integer>> permuteNodes(ArrayList<ArrayList<Integer>> nodes)
        {
            ArrayList<ArrayList<ArrayList<Integer>>> permutationsForChildren = new ArrayList<>();

            ArrayList<Integer> currentPermutation = new ArrayList<>();
            ArrayList<ArrayList<Integer>> finalPermutations = new ArrayList<>();

            // Build List with the permutation of their children, for each Node
            for (ArrayList<Integer> children : nodes)
            {
                permutationsForChildren.add(permuteList(children));
            }

            permuteNodes(permutationsForChildren, currentPermutation, finalPermutations, 0);

            return finalPermutations;
        }

        // We want to permute all these list so...
        // a[0] + a1[0] + a2[0] + a3[0]
        // a[0] + a1[0] + a2[0] + a3[1]
        // and so on...
        // The pm structure is structured like so
        // Parent Child -> i'th Permutation Reference -> i'th Permutation
        // We will be returning all possible permuted lists, as shown above.
        private void permuteNodes(ArrayList<ArrayList<ArrayList<Integer>>> permutationsForChildren, 
                                                ArrayList<Integer> currentPermutation, 
                                                ArrayList<ArrayList<Integer>> finalPermutations, 
                                                int node)
        {
            // Stop if we permuted for each Parent node
            if (node == permutationsForChildren.size())
            {
                finalPermutations.add(new ArrayList<>(currentPermutation));
                return;
            }

            // Get all permutations for current Node
            for (ArrayList<Integer> permutation : permutationsForChildren.get(node))
            {
                ArrayList<Integer> currPermCopy = new ArrayList<>(currentPermutation);
                currPermCopy.addAll(permutation);

                permuteNodes(permutationsForChildren, currPermCopy, finalPermutations, node + 1);
            }
        }

        private ArrayList<ArrayList<Integer>> permuteList(ArrayList<Integer> list)
        {
        	ArrayList<ArrayList<Integer>> permuted = new ArrayList<>();
        	permuteList(permuted, list, list.size());

            return permuted;
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

        public void dfs(int node)
        {
            boolean [] visited = new boolean[n];

            ArrayList<Integer> path = new ArrayList<>(n);

            visited[node] = true;
            path.add(0);

            System.out.println("DFS Traversal...");
            dfs(visited, path, node);
        }

        private boolean dfs(boolean [] visited, ArrayList<Integer> path, int node)
        {
            if (path.size() == n)
            {
                System.out.println("\tFound path...");
                System.out.print("\t");
                printPath(path);
                System.out.println();
                return true;
            }

            for (int i = 0; i < n; i++)
            {
                if (matrix[node][i] && !visited[i])
                {
                    visited[i] = true;
                    path.add(i);

                    if (dfs(visited, path, i))
                    {
                        return true;
                    }

                    visited[i] = false;
                }
            }

            return false;
        }

        // This method will print out all possible dfs traversals
        // 
        public void dfsAll(int startNode)
        {
            boolean [] visited = new boolean[n];
            int [] path = new int[n];
            Arrays.fill(path, -1);

            Stack<Integer> s = new Stack<Integer>();

            visited[startNode] = true;
            s.push(startNode);


            System.out.println("All DFS Traversal Iteratively...");

            // To traverse, we must keep adding to the stack
            boolean flag = false;
            while (!s.isEmpty())
            {
                flag = false;
                int node = s.peek();
                printNode(node);
                System.out.println();
                for(ArrayList<Integer> permutation : permuteList(list.get(node)))
                {
                    for (Integer aNode : permutation)
                    {
                        if (!visited[aNode])
                        {
                            visited[aNode] = true;
                            s.push(aNode);
                            flag = true;
                        }
                    }
                }

                if (!flag)
                {
                    s.pop();
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

	// We want to have DFS iteratively with all traversals ** HARD **
	
	// Check if graph is bipartite
        // Can simply adding to a different array (have two arrays, insert at opposite arrays at each iteration),
        // work? In the case of an empty graph, it is necessarily true? The only time that a graph is NOT bipartite
        // is when there are no ways to construct two Sets of Nodes where they are not adjacent to each other.
        //
        // It seems that BFS must work. For every iteration of BFS, we are visiting Nodes that are part of the
        // Set opposing our current Node. If for every iteration, we do BFS traversal and add the visited Nodes
        // to the corresponding Set, and none of those Nodes are adjacent each other, we reach our solution.
        //
        // Important to note:
        // It may be helpful to know that we may only return false if we notice that two Nodes are adjacent in a BFS traversal.
        
        public boolean isBipartite(int startNode)
        {
            boolean [] used = new boolean[n];
            // we need to keep track of what Nodes are at either side of the graph
            boolean [][] partition = new boolean[n][n];

            Queue<Integer> q = new ArrayDeque<>();

            used[startNode] = true;
            partition[0][startNode] = true;
            q.add(startNode);

            int part = 0;
            int itr = 0;
            while (itr < n)
            {
                int qSize = q.size();
                int [] nodes = new int[qSize];

                for (int i = 0; i < qSize; i++)
                {
                    itr++;
                    nodes[i] = q.remove();
                    partition[part][nodes[i]] = true;
                }

                // Check if any other Node in the partition is adjacent to this
                for (int j = 0; j < qSize; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        // Find any element that have been inserted
                        if (partition[part][i] && i != nodes[j])
                        {
                            // Check if there is adjacency
                            if (matrix[nodes[j]][i] || matrix[i][nodes[j]])
                            {
                                return false;
                            }
                        }
                    }
                }

                for (int j = 0; j < qSize; j++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        if (matrix[nodes[j]][i] && !used[i])
                        {
                            q.add(i);
                            used[i] = true;
                        }
                    }
                }

                part = (part + 1) % 2;
            }

            return true;
        }

	
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
            if (path[i] == -1) return;
			System.out.print(value[path[i]] + " -> ");
		}
		
		if (path.length > 1)
		{
			System.out.println(value[path[path.length-1]]);
		}
	}

    private void printPath(ArrayList<Integer> path)
    {
        for (int i = 0; i < path.size() - 1; i++)
        {
            System.out.print(value[path.get(i)] + " -> ");
        }
        
        if (path.size() > 1)
        {
            System.out.println(value[path.get(path.size()-1)]);
        }
    }

    private void printPathAsSet(ArrayList<Integer> path)
    {
        System.out.print("[");
        for (int i = 0; i < path.size() - 1; i++)
        {
            System.out.print(value[path.get(i)] + ", ");
        }
        
        if (path.size() >= 1)
        {
            System.out.print(value[path.get(path.size()-1)]);
        }

        System.out.println("]");
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
            g.topologicalSortIteratively(0);
            g.bfsIteratively(0);
            //g.bfs(0);
            g.dfs(0);
            g.dfsAll(0);
            System.out.println("isBipartite() " + g.isBipartite(0));
    }
}
