package chromo_activities;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import chromosome.InitializeChromosome;
import clustered.Cluster;
import graph.DijkstraHeap;
import graph.GraphMethods;
import iofile.ReadWriteFile;

public class Crossover {
	GraphMethods graph = new GraphMethods();

	private double[][] findSpanningTreeBetweenClusters(double[][] par_1, int num_Vertex, int num_Cluster,
			int[][] vertex_In_Cluster) {
		double[][] T = new double[num_Cluster][num_Cluster];
		// Khoi tao
		for (int i = 0; i < num_Cluster; i++) {
			for (int j = i; j < num_Cluster; j++) {
				T[i][j] = 0;
				T[j][i] = 0;
			}
		}

		int num_Vertex_in_Cluster_1 = 0;
		int num_Vertex_in_Cluster_2 = 0;
		for (int i = 0; i < num_Cluster; i++) {
			num_Vertex_in_Cluster_1 = vertex_In_Cluster[i].length;
			for (int j = 0; j < num_Vertex_in_Cluster_1; j++) {
				for (int k = i + 1; k < num_Cluster; k++) {
					num_Vertex_in_Cluster_2 = vertex_In_Cluster[k].length;
					for (int t = 0; t < num_Vertex_in_Cluster_2; t++) {
						if (par_1[vertex_In_Cluster[i][j]][vertex_In_Cluster[k][t]] > 0) {
							T[i][k] = par_1[vertex_In_Cluster[i][j]][vertex_In_Cluster[k][t]];
							T[k][i] = T[i][k];
						}
					}
				}
			}
		}
		return T;
	}

	/************************ Lai ghep bang Prim *********************/

	private double[][] getUnionOfClusters(double[][] par1, double[][] par2, int num_Vertex, Random rnd) {
		double[][] G_cr = new double[num_Vertex][num_Vertex];
		InitializeChromosome init_Chrome = new InitializeChromosome();
		// 01. Tao do thi Gcr
		for (int i = 0; i < num_Vertex; i++) {
			for (int j = 0; j < num_Vertex; j++) {
				G_cr[i][j] = par1[i][j] + par2[i][j];
				if (G_cr[i][j] > 1)
					G_cr[i][j] = 1;
			}
		}
		return init_Chrome.primRST(num_Vertex, G_cr, rnd);
	}

	public double[][] crossoverByPrimRST(double[][] par_1, double[][] par_2, int num_Vertex, int num_Cluster,
			int[][] vertex_In_Cluster, Random rnd) {
		double[][] T = new double[num_Vertex][num_Vertex];
		InitializeChromosome init_Chrome = new InitializeChromosome();

		// Khoi tao
		for (int i = 0; i < num_Vertex; i++) {
			for (int j = i; j < num_Vertex; j++) {
				T[i][j] = 0;
				T[j][i] = 0;
			}
		}

		// Ap dung PrimRST crossover cho tung cluster
		double[][] weight_Cluster_Matrix_Par_1, weight_Cluster_Matrix_Par_2;
		double[][] spanning_Tree_of_Cluster;
		int num_Vertex_in_Cluster = 0;
		for (int i = 0; i < num_Cluster; i++) {
			num_Vertex_in_Cluster = vertex_In_Cluster[i].length;
			weight_Cluster_Matrix_Par_1 = init_Chrome.createWeightMatrixForCluster(num_Vertex, num_Vertex_in_Cluster,
					par_1, vertex_In_Cluster[i]);
			weight_Cluster_Matrix_Par_2 = init_Chrome.createWeightMatrixForCluster(num_Vertex, num_Vertex_in_Cluster,
					par_2, vertex_In_Cluster[i]);
			spanning_Tree_of_Cluster = getUnionOfClusters(weight_Cluster_Matrix_Par_1, weight_Cluster_Matrix_Par_2,
					num_Vertex_in_Cluster, rnd);
			// Chuyen ra cay khung cua do thi G
			for (int k = 0; k < num_Vertex_in_Cluster; k++) {
				for (int j = 0; j < num_Vertex_in_Cluster; j++) {
					if (spanning_Tree_of_Cluster[k][j] > 0) {
						T[vertex_In_Cluster[i][k]][vertex_In_Cluster[i][j]] = spanning_Tree_of_Cluster[k][j];
					}
				}
			}
		}
		// Tạo cây khung cho đại diện các nhóm
		// 01. Tim cạnh liên kết các nhóm của cha mẹ => các cạnh của cha mẹ có
		// thể liên kết tới các đỉnh khác nhau của cùng 1 nhóm =>
		weight_Cluster_Matrix_Par_1 = findSpanningTreeBetweenClusters(par_1, num_Vertex, num_Cluster,
				vertex_In_Cluster);
		weight_Cluster_Matrix_Par_2 = findSpanningTreeBetweenClusters(par_2, num_Vertex, num_Cluster,
				vertex_In_Cluster);

		// 02. Tạo các thể con cho các cạnh nay
		spanning_Tree_of_Cluster = getUnionOfClusters(weight_Cluster_Matrix_Par_1, weight_Cluster_Matrix_Par_2,
				num_Cluster, rnd);
		// Chuyen ra cay khung cua do thi G
		int[] idx_Cluster = new int[num_Cluster];
		for (int i = 0; i < num_Cluster; i++) {
			int vt = rnd.nextInt(vertex_In_Cluster[i].length);
			idx_Cluster[i] = vertex_In_Cluster[i][vt];
		}
		for (int k = 0; k < num_Cluster; k++) {
			for (int j = 0; j < num_Cluster; j++) {
				if (spanning_Tree_of_Cluster[k][j] > 0) {
					T[idx_Cluster[k]][idx_Cluster[j]] = spanning_Tree_of_Cluster[k][j];
				}
			}
		}
		return T;
	}

	/************************** End region **********************************/

	/******************** Lai ghep bang BFS ****************************/

	public double[][] BFSCrossover(double[][] father, double[][] mother, int num_vertex, int startVertex) {

		double[][] parentGroup = new double[num_vertex][num_vertex];
		double[][] spanningTree = new double[num_vertex][num_vertex];
		// addition two matrix of father and mother

		for (int i = 0; i < num_vertex; i++) {
			for (int j = 0; j < num_vertex; j++) {
				parentGroup[i][j] = father[i][j] + mother[i][j];
			}
		}
		spanningTree = graph.BFSTree(parentGroup, num_vertex, startVertex);

		return spanningTree;
	}

	/**
	 * 
	 * @param father
	 * @param mother
	 * @param num_vertex
	 * @param clusters
	 * @return
	 */
	public double[][] crossoverByBFS(double[][] father, double[][] mother, int num_vertex, int[][] vertexInCluster,
			Random rnd) {

		double[][] Tree = new double[num_vertex][num_vertex];
		double[][] clusterWeightMatrix1, clusterWeightMatrix2;
		double[][] spanningTreeOfCluster;
		int numberClusterVertex = 0;
		int numberOfCluster = vertexInCluster.length;
		// each Cluster is presented by one vertex in that cluster
		int[] presentationVertex = new int[numberOfCluster];
		int[] indexPresentationVertex = new int[numberOfCluster];
		// choose randomly a vertex in each cluster

		for (int i = 0; i < numberOfCluster; i++) {
			numberClusterVertex = vertexInCluster[i].length;
			// we choose root is presentation vertex for cluster which contains
			// root vertex
			if (i == this.idCluster()) {
				for (int j = 0; j < vertexInCluster[i].length; j++) {
					if (vertexInCluster[i][j] == ReadWriteFile.root) {
						indexPresentationVertex[i] = j;
					}
				}
			}

			else {
				indexPresentationVertex[i] = rnd.nextInt(numberClusterVertex);
			}
			presentationVertex[i] = (int) vertexInCluster[i][indexPresentationVertex[i]];

			InitializeChromosome initializeChromosome = new InitializeChromosome();
			clusterWeightMatrix1 = initializeChromosome.createWeightMatrixForCluster(num_vertex,
					vertexInCluster[i].length, father, vertexInCluster[i]);
			clusterWeightMatrix2 = initializeChromosome.createWeightMatrixForCluster(num_vertex,
					vertexInCluster[i].length, mother, vertexInCluster[i]);
			spanningTreeOfCluster = BFSCrossover(clusterWeightMatrix1, clusterWeightMatrix2, numberClusterVertex,
					indexPresentationVertex[i]);

			// copy to graph
			for (int j = 0; j < numberClusterVertex; j++) {
				for (int k = 0; k < numberClusterVertex; k++) {
					Tree[vertexInCluster[i][j]][vertexInCluster[i][k]] = spanningTreeOfCluster[j][k];
				}
			}

		}
		// build the vertex representation for each cluster
		clusterWeightMatrix1 = findSpanningTreeBetweenClusters(father, num_vertex, numberOfCluster, vertexInCluster);
		clusterWeightMatrix2 = findSpanningTreeBetweenClusters(mother, num_vertex, numberOfCluster, vertexInCluster);

		spanningTreeOfCluster = getUnionOfClusters(clusterWeightMatrix1, clusterWeightMatrix2, numberOfCluster, rnd);

		// convert to spanning tree of Graph
		int[] indexCluster = new int[numberOfCluster];
		for (int i = 0; i < numberOfCluster; i++) {
			int position = indexPresentationVertex[i];
			indexCluster[i] = vertexInCluster[i][position];
			for (int j = 0; j < numberOfCluster; j++) {
				for (int k = 0; k < numberOfCluster; k++) {
					if (spanningTreeOfCluster[j][k] > 0) {
						Tree[indexCluster[j]][indexCluster[k]] = spanningTreeOfCluster[j][k];
					}
				}
			}
		}
		return Tree;
	}

	public int idCluster() {
		int id = 0;
		boolean b = true;
		for (int i = 0; i < ReadWriteFile.clusters.size(); i++) {
			int clusterVertex = ReadWriteFile.clusters.get(i).getCluster().size();
			for (int j = 0; j < clusterVertex; j++) {
				if (ReadWriteFile.root == ReadWriteFile.clusters.get(i).getCluster().get(j)) {
					id = i;
					b = false;
					break;
				}
			}
			if (b == false)
				break;
		}
		return id;
	}

	/**************** End region ************************/

	/*************** Lai ghep bang Dijkstra *******************/
	
	public double[][] crossOverByDijkstra(double[][] par1, double[][] par2, int[][] vertexInCluster, Random rnd) {
		DijkstraHeap dijkstra = new DijkstraHeap();
		ArrayList<Cluster> p1 = dijkstra.convertArrayToCluster(par1, vertexInCluster);
		ArrayList<Cluster> p2 = dijkstra.convertArrayToCluster(par2, vertexInCluster);
		ArrayList<Cluster> child = new ArrayList<Cluster>();
		for (int i = 0; i < p1.size(); i++) {
			int start = vertexInCluster[i][rnd.nextInt(vertexInCluster[i].length)];
			child.add(dijkstra.Al_Dijkstra_Heap(p1.get(i), p2.get(i), start));
		}
		double[][] T = new double[ReadWriteFile.numberOfCity][ReadWriteFile.numberOfCity];
		for (int i = 0; i < ReadWriteFile.numberOfCity; i++) {
			Arrays.fill(T[i], 0);
		}
		for (int i = 0; i < child.size(); i++) {
			T = dijkstra.convertClusterToArray(child.get(i), T);
		}
		for (int i = 0; i < vertexInCluster.length; i++) {
			for (int j = i + 1; j < vertexInCluster.length; j++) {
				for (int k = 0; k < vertexInCluster[i].length; k++) {
					for (int l = 0; l < vertexInCluster[j].length; l++) {
						if (par1[vertexInCluster[i][k]][vertexInCluster[j][l]] == 1) {
							T[vertexInCluster[i][k]][vertexInCluster[j][l]] = 1;
							T[vertexInCluster[j][l]][vertexInCluster[i][k]] = 1;
						}
					}
				}
			}
		}
		return T;
	}

}
