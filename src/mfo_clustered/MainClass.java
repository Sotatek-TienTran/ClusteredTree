package mfo_clustered;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class MainClass {
	ReadWriteFile ioFile = new ReadWriteFile();
	InitializeChromosome initChromo = new InitializeChromosome();
	Evaluate eva = new Evaluate();
	Chromosome chromo = new Chromosome();
	Crossover cross = new Crossover();
	Mutation mutation = new Mutation();
	private static double crossOverRate = 0.9;
	private static double mutationRate = 0.9;
	static Random rnd = new Random(1);
	public int maxGroup;
	public static int defaultPopLength = 100;
	ArrayList<Chromosome> population = new ArrayList<Chromosome>();

	public static void main(String[] args) {
		MainClass mainClass = new MainClass();
		int numOfCity = ReadWriteFile.numberOfCity;
		int numOfCluster = ReadWriteFile.numberOfCluster;
		int vertexInCluster[][] = ReadWriteFile.vertexInCluster;
		ArrayList<Chromosome> pop = init();

		System.out.print(pop.get(0).getFactorialCost()[1] + " ");
		sortByCostIndex(pop, 0);
		System.out.println(pop.get(0).getFactorialCost()[0]);
		for (int i = 0; i < 50; i++) {
			ArrayList<Chromosome> tempPop = new ArrayList<Chromosome>();
			while (tempPop.size() < defaultPopLength) {
				double r = 0 + (1 - 0) * rnd.nextDouble();
				int par1 = rnd.nextInt(defaultPopLength);
				int par2 = rnd.nextInt(defaultPopLength);
				if (pop.get(par1).getSkillFactor() == pop.get(par2).getSkillFactor() || r < crossOverRate) {
					Chromosome child = new Chromosome();
					child.setEdgesMatrix(mainClass.cross.primRSTClusterCrossover(pop.get(par1).getEdgesMatrix(),
							pop.get(par1).getEdgesMatrix(), numOfCity, numOfCluster, vertexInCluster, rnd));
					tempPop.add(child);
				} else {
					pop.get(par1)
							.setEdgesMatrix(mainClass.mutation.edgeClusteredTreeMutation(pop.get(par1).getEdgesMatrix(),
									numOfCity, mutationRate, numOfCluster, vertexInCluster, rnd));
					tempPop.add(pop.get(par1));
					pop.get(par2)
							.setEdgesMatrix(mainClass.mutation.edgeClusteredTreeMutation(pop.get(par2).getEdgesMatrix(),
									numOfCity, mutationRate, numOfCluster, vertexInCluster, rnd));
					tempPop.add(pop.get(par2));
				}
			}
			if (tempPop.size() > defaultPopLength) {
				tempPop.remove(rnd.nextInt(defaultPopLength));
			}
			calculate(tempPop);
			evaluatePopulation(tempPop);
			//drawResult(tempPop.get(0).getEdgesMatrix());
			ArrayList<Chromosome> newPop = new ArrayList<Chromosome>();
			for (int j = 0; j < defaultPopLength / 2; j++) {
				newPop.add(pop.get(j));
				newPop.add(tempPop.get(j));
			}
			pop = newPop;
			
		}
		calculate(pop);
		evaluatePopulation(pop);
		System.out.print(pop.get(0).getFactorialCost()[1] + " ");
		// drawResult(mainClass.eva.decodingMFOVertexInSubGraph(pop.get(0).getEdgesMatrix(),
		// 14, 14));
		sortByCostIndex(pop, 0);
		System.out.println(pop.get(0).getFactorialCost()[0]);
		drawResult(pop.get(0).getEdgesMatrix());
//		Crossover c = new Crossover();
//		double temp[][] = new double[numOfCity][numOfCity];
//		temp = c.ClusterBFSCrossover(pop.get(0).getEdgesMatrix(), pop.get(1).getEdgesMatrix(), numOfCity, vertexInCluster, rnd);
//		printArray(temp);
	}

	public static void evaluatePopulation(ArrayList<Chromosome> pop) {
		sortByCostIndex(pop, 0);
		for (int i = 0; i < defaultPopLength; i++) {
			pop.get(i).setFactorialRank(i + 1, 0);
			// System.out.println(pop.get(i).getFactorialCost()[0] + " ");
		}
		sortByCostIndex(pop, 1);
		for (int i = 0; i < defaultPopLength; i++) {
			pop.get(i).setFactorialRank(i + 1, 1);
		}
		for (int i = 0; i < defaultPopLength; i++) {
			if (pop.get(i).getFactorialRank()[0] < pop.get(i).getFactorialRank()[1]) {
				pop.get(i).setSkillFactor(0);
				pop.get(i).setScalarFitness(1.0 / (pop.get(i).getFactorialRank()[0]));
			} else {
				pop.get(i).setSkillFactor(1);
				pop.get(i).setScalarFitness(1.0 / (pop.get(i).getFactorialRank()[1]));
			}

		}
		Collections.sort(pop, ChromosomeCmp.compareByScalarFitness);
	}

	public static void drawResult(double[][] weightMatrix) {
		City cities[] = ReadWriteFile.cities;
		for (int i = 0; i < weightMatrix.length; i++) {
			
			for (int j = 0; j < weightMatrix[i].length; j++) {
				if (weightMatrix[i][j] > 0) {
					DrawLines d = new DrawLines();
					d.draw(cities[i].getX(), cities[i].getY(), cities[j].getX(), cities[j].getY(), i + 1 + "",
							j + 1 + "");
				}
			}
		}
	}

	public static ArrayList<Chromosome> init() {
		MainClass mc = new MainClass();
		int numOfCity = ReadWriteFile.numberOfCity;
		int numOfCluster = ReadWriteFile.numberOfCluster;
		double weightMatrix[][] = ReadWriteFile.weightMatrix;
		int vertexInCluster[][] = ReadWriteFile.vertexInCluster;
		int maxGroupValue = 0;
		for (int i = 0; i < numOfCluster; i++) {
			if (maxGroupValue < vertexInCluster[i].length) {
				maxGroupValue = vertexInCluster[i].length;
			}
		}
		mc.maxGroup = maxGroupValue;
		for (int i = 0; i < defaultPopLength; i++) {
			Chromosome c = new Chromosome();
			c.setEdgesMatrix(
					mc.initChromo.primRSTForClusteredTree(numOfCity, numOfCluster, weightMatrix, vertexInCluster, rnd));

			mc.population.add(c);
		}

		calculate(mc.population);
		return mc.population;
	}

	public static void calculate(ArrayList<Chromosome> pop) {
		MainClass mc = new MainClass();
		int numOfCity = ReadWriteFile.numberOfCity;
		int numOfCluster = ReadWriteFile.numberOfCluster;
		double weightMatrix[][] = ReadWriteFile.weightMatrix;
		int vertexInCluster[][] = ReadWriteFile.vertexInCluster;
		int startVertex = ReadWriteFile.root;
		int maxGroupValue = 0;
		for (int i = 0; i < numOfCluster; i++) {
			if (maxGroupValue < vertexInCluster[i].length) {
				maxGroupValue = vertexInCluster[i].length;
			}
		}
		for (int i = 0; i < defaultPopLength; i++) {
			double[] temp = new double[2];
			temp[0] = mc.eva.evaluation(pop.get(i).getEdgesMatrix(), weightMatrix, numOfCity, startVertex);
			temp[1] = mc.eva.evaluation(
					mc.eva.decodingMFOVertexInSubGraph(pop.get(i).getEdgesMatrix(), maxGroupValue, maxGroupValue),
					weightMatrix, maxGroupValue, rnd.nextInt(maxGroupValue));
			pop.get(i).setFactorialCost(temp);
		}

	}

	public static void sortByCostIndex(ArrayList<Chromosome> pop, int index) {
		for (int i = 0; i < defaultPopLength; i++) {
			pop.get(i).cost = pop.get(i).getFactorialCost()[index];
		}
		Collections.sort(pop, ChromosomeCmp.compareByFactorialCost);
	}

	public static void printArray(double[][] arr) {
		for (int i = 0; i < arr.length; i++) {
			for (int j = 0; j < arr[i].length; j++) {
				System.out.print(arr[i][j] + " ");
			}
			System.out.println();
		}

	}
}
