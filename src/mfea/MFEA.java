package mfea;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import chromo_activities.Crossover;
import chromo_activities.Decode;
import chromo_activities.Mutation;
import chromosome.Chromosome;
import chromosome.InitializeChromosome;
import evaluate.Evaluate;
import graph.GraphMethods;
import iofile.ReadWriteFile;
import chromosome.ChromosomeCmp;
import clustered.City;
import draw.DrawLines;

public class MFEA {
	ReadWriteFile ioFile = new ReadWriteFile();
	InitializeChromosome initChromo = new InitializeChromosome();
	Evaluate eva = new Evaluate();
	Chromosome chromo = new Chromosome();
	Crossover cross = new Crossover();
	Mutation mutation = new Mutation();
	GraphMethods graph = new GraphMethods();
	Decode decode = new Decode();
	static double crossOverRate = 0.5;
	static double mutationRate = 0.05;
	static Random rnd;
	public int maxGroup;
	public static int defaultPopLength = 100;
	public static int numberOfGeneration = 50000;
	ArrayList<Chromosome> population = new ArrayList<Chromosome>();

	public static void main(String[] args) {
		ArrayList<Chromosome> pop = new ArrayList<Chromosome>();
		// rnd = new Random(1);
		for (int i = 0; i < 30; i++) {
			rnd = new Random(i);
//			 ReadWriteFile.clusterReadFiles("Instance/10eil51.clt");
			ReadWriteFile.clusterReadFiles("Type_1_Small/" + args[0] + ".clt");

			PrintWriter pw1 = null;
			// rnd = new Random(i);
			// pop = GA.geneticAlgorithm(numberOfGeneration, rnd);
			pop = MFO(numberOfGeneration, rnd);
			// System.out.println(pop.get(0).getFactorialCost()[0]);
			// result = result + pop.get(0).getFactorialCost()[0];
			try {
				pw1 = new PrintWriter(new FileWriter(new File(args[1]), true));
			} catch (IOException e) {
				e.printStackTrace();
			}
			pw1.print(pop.get(0).getFactorialCost()[0]+"\t");
			sortByCostIndex(pop, 1);
			pw1.println(pop.get(0).getFactorialCost()[1]);
			pw1.close();
//			for (int j = 0; j < defaultPopLength; j++) {
//				System.out.print(pop.get(j).getFactorialCost()[0] + " ");
//			}

		}
	}

	public static ArrayList<Chromosome> MFO(int numberOfGeneration, Random rnd) {
		MFEA mFEA = new MFEA();
		int numOfCity = ReadWriteFile.numberOfCity;
		int numOfCluster = ReadWriteFile.numberOfCluster;
		int vertexInCluster[][] = ReadWriteFile.vertexInCluster;
		ArrayList<Chromosome> pop = init();
		for (int i = 0; i < numberOfGeneration; i++) {
			ArrayList<Chromosome> tempPop = new ArrayList<Chromosome>();
			while (tempPop.size() < defaultPopLength) {
				double r = 0 + (1 - 0) * rnd.nextDouble();
				int par1 = rnd.nextInt(defaultPopLength);
				int par2;
				do {
					par2 = rnd.nextInt(defaultPopLength);
				} while (par2 == par1);
				Chromosome temp1 = new Chromosome();
				Chromosome temp2 = new Chromosome();
				temp1.setEdgesMatrix(pop.get(par1).getEdgesMatrix());
				temp2.setEdgesMatrix(pop.get(par2).getEdgesMatrix());
				if ((pop.get(par1).getSkillFactor() == pop.get(par1).getSkillFactor() || r < crossOverRate)) {
					Chromosome child = new Chromosome();
					// child.setEdgesMatrix(mFEA.cross.crossoverByPrimRST(temp1.getEdgesMatrix(),
					// temp2.getEdgesMatrix(),
					// numOfCity, numOfCluster, vertexInCluster, rnd));
//					child.setEdgesMatrix(mFEA.cross.crossoverByBFS(temp1.getEdgesMatrix(), temp2.getEdgesMatrix(),
//							numOfCity, vertexInCluster, rnd));
					child.setEdgesMatrix(mFEA.cross.crossOverByDijkstra(temp1.getEdgesMatrix(), temp2.getEdgesMatrix(), vertexInCluster, rnd));
					child.setEdgesMatrix(mFEA.mutation.mutationTree(child.getEdgesMatrix(), numOfCity, mutationRate,
							numOfCluster, vertexInCluster, rnd));
					tempPop.add(child);
				} else {
					temp1.setEdgesMatrix(mFEA.mutation.mutationTree(temp1.getEdgesMatrix(), numOfCity, mutationRate,
							numOfCluster, vertexInCluster, rnd));
					tempPop.add(temp1);
					temp2.setEdgesMatrix(mFEA.mutation.mutationTree(temp2.getEdgesMatrix(), numOfCity, mutationRate,
							numOfCluster, vertexInCluster, rnd));
					tempPop.add(temp2);
				}
				if (pop.get(par1).equals(temp1)) {
					System.out.println(0);
				}
			}
			if (tempPop.size() > defaultPopLength) {
				tempPop.remove(rnd.nextInt(tempPop.size()));
			}
			calculate(tempPop);
			evaluatePopulation(tempPop);
			//
			// sortByCostIndex(tempPop, 0);
			// tempPop.remove(rnd.nextInt(tempPop.size()));
			ArrayList<Chromosome> newPop = new ArrayList<Chromosome>();
			for (int j = 0; j < defaultPopLength - 5; j++) {
				newPop.add(tempPop.get(j));
			}
			for (int j = 0; j < 5; j++) {
				newPop.add(pop.get(j));
			}
			// newPop.add(pop.get(0));
			pop.clear();
			for (int j = 0; j < defaultPopLength; j++) {
				pop.add(newPop.get(j));
			}
			calculate(pop);
			evaluatePopulation(pop);
//			if (i % 100 == 0) {
//				System.out.println(pop.get(0).getFactorialCost()[0] + " " + pop.get(1).getFactorialCost()[0]);
//			}
		}
		return pop;
	}

	public static ArrayList<Chromosome> init() {
		MFEA mc = new MFEA();
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
		evaluatePopulation(mc.population);
		return mc.population;
	}

	public static void calculate(ArrayList<Chromosome> pop) {
		MFEA mc = new MFEA();
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
					mc.decode.decodingMFOVertexInSubGraph(pop.get(i).getEdgesMatrix(), maxGroupValue, maxGroupValue),
					weightMatrix, maxGroupValue, rnd.nextInt(maxGroupValue));
			pop.get(i).setFactorialCost(temp);
		}
	}

	public static void evaluatePopulation(ArrayList<Chromosome> pop) {
		sortByCostIndex(pop, 1);
		for (int i = 0; i < defaultPopLength; i++) {
			pop.get(i).setFactorialRank(i + 1, 1);
			// System.out.println(pop.get(i).getFactorialCost()[1]);
		}
		sortByCostIndex(pop, 0);
		for (int i = 0; i < defaultPopLength; i++) {
			pop.get(i).setFactorialRank(i + 1, 0);
			// System.out.println(pop.get(i).getFactorialCost()[0] + " ");
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
		// Collections.sort(pop, ChromosomeCmp.compareByScalarFitness);
	}

	public static void sortByCostIndex(ArrayList<Chromosome> pop, int index) {
		for (int i = 0; i < defaultPopLength; i++) {
			pop.get(i).setCost(pop.get(i).getFactorialCost()[index]);
		}
		Collections.sort(pop, ChromosomeCmp.compareByFactorialCost);
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

}
