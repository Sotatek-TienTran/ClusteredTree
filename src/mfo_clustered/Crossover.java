package mfo_clustered;

import java.util.Random;

public class Crossover {
	InitializeChromosome init_Chrom = new InitializeChromosome();

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

	public double[][] primRSTClusterCrossover(double[][] par_1, double[][] par_2, int num_Vertex, int num_Cluster, int[][] vertex_In_Cluster, Random rnd)
     {
         double[][] T = new double[num_Vertex][num_Vertex];
         InitializeChromosome init_Chrome = new InitializeChromosome();

         //Khoi tao
         for (int i = 0; i < num_Vertex; i++)
         {
             for (int j = i; j < num_Vertex; j++)
             {
                 T[i][j] = 0;
                 T[j][i] = 0;
             }
         }


         //Ap dung PrimRST crossover cho tung cluster
         double[][] weight_Cluster_Matrix_Par_1, weight_Cluster_Matrix_Par_2;
         double[][] spanning_Tree_of_Cluster;
         int num_Vertex_in_Cluster = 0;
         for(int i = 0; i < num_Cluster; i ++)
         {
             num_Vertex_in_Cluster = vertex_In_Cluster[i].length;
             weight_Cluster_Matrix_Par_1 = init_Chrome.createWeightMatrixForCluster(num_Vertex, num_Vertex_in_Cluster, par_1, vertex_In_Cluster[i]);
             weight_Cluster_Matrix_Par_2 = init_Chrome.createWeightMatrixForCluster(num_Vertex, num_Vertex_in_Cluster, par_2, vertex_In_Cluster[i]);
             spanning_Tree_of_Cluster = primRSTCrossover(weight_Cluster_Matrix_Par_1, weight_Cluster_Matrix_Par_2, num_Vertex_in_Cluster, rnd);
             //Chuyen ra cay khung cua do thi G
             for (int k = 0; k < num_Vertex_in_Cluster; k++)
             {
                 for (int j = 0; j < num_Vertex_in_Cluster; j++)
                 {
                     if (spanning_Tree_of_Cluster[k][j] > 0)
                     {
                         T[vertex_In_Cluster[i][k]][vertex_In_Cluster[i][j]] = spanning_Tree_of_Cluster[k][j];
                     }
                 }
             }
         }
         //Tạo cây khung cho đại diện các nhóm 
         //01. Tim cạnh liên kết các nhóm của cha mẹ => các cạnh của cha mẹ có thể liên kết tới các đỉnh khác nhau của cùng 1 nhóm =>  
         weight_Cluster_Matrix_Par_1 = findSpanningTreeBetweenClusters(par_1, num_Vertex, num_Cluster, vertex_In_Cluster);
         weight_Cluster_Matrix_Par_2 = findSpanningTreeBetweenClusters(par_2, num_Vertex, num_Cluster, vertex_In_Cluster);
         
         //02. Tạo các thể con cho các cạnh nay
         spanning_Tree_of_Cluster = primRSTCrossover(weight_Cluster_Matrix_Par_1, weight_Cluster_Matrix_Par_2, num_Cluster, rnd);
         //Chuyen ra cay khung cua do thi G
         int[] idx_Cluster = new int[num_Cluster];
         for (int i = 0; i < num_Cluster; i++)
         {
             int vt = rnd.nextInt(vertex_In_Cluster[i].length);
             idx_Cluster[i] = vertex_In_Cluster[i][vt];
         }
         for (int k = 0; k < num_Cluster; k++)
         {
             for (int j = 0; j < num_Cluster; j++)
             {
                 if (spanning_Tree_of_Cluster[k][j] > 0)
                 {
                     T[idx_Cluster[k]][idx_Cluster[j]] = spanning_Tree_of_Cluster[k][j];
                 }
             }
         }
         return T;
     }

	private double[][] primRSTCrossover(double[][] par1, double[][] par2, int num_Vertex, Random rnd) {
		double[][] G_cr = new double[num_Vertex][num_Vertex];
		InitializeChromosome init_Chrome = new InitializeChromosome();
		// 01. Tao do thi Gcr
		for (int i = 0; i < num_Vertex; i++) {
			for (int j = 0; j < num_Vertex; j++) {
				G_cr[i][j] = par1[i][j] + par2[i][j];
			}
		}
		return init_Chrome.primRST(num_Vertex, G_cr, rnd);
	}

	public Chromosome MFOCrossover(double crss_Para, Chromosome par_1, Chromosome par_2, int num_Genens, Random rnd){
		Chromosome child = new Chromosome();
		child.setEdgesMatrix(primRSTCrossover(par_1.getEdgesMatrix(), par_2.getEdgesMatrix(), num_Genens, rnd));
		return child;
	}
}
