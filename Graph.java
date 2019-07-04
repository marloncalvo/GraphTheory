import java.util.ArrayList;
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
	
	// We want to have BFS iteratively
	
	// We want to have DFS for graph recursively
	
	// We want to have DFS iteratively
	
	// Do the same for all, with Lists
	
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
	}
}
