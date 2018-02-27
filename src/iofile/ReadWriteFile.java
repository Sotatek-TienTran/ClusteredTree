package iofile;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import clustered.Cluster;
import clustered.City;

public class ReadWriteFile {
	public static int numberOfCity;
	public static double weightMatrix[][];
	public static ArrayList<Cluster> clusters;
	public static int numberOfCluster;
	public static int root;
	public static int[][] vertexInCluster;
	public static City[] cities;
	public static int maxGroupIndex;

	public static void clusterReadFiles(String fileName) {
		clusters = new ArrayList<Cluster>();
		BufferedReader br = null;
		try {
			String sCurrentLine = null;
			br = new BufferedReader(new FileReader(fileName));
			for (int j = 0; j < 3; j++) {
				sCurrentLine = br.readLine();
			}
			String[] str = sCurrentLine.split(": ");
			numberOfCity = Integer.parseInt(str[1]);

			weightMatrix = new double[numberOfCity][numberOfCity];
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			numberOfCluster = Integer.parseInt(str[1]);
			sCurrentLine = br.readLine();
			sCurrentLine = br.readLine();

			// read the detail of the City
			cities = new City[numberOfCity];
			for (int j = 0; j < numberOfCity; j++) {
				sCurrentLine = br.readLine();
				str = sCurrentLine.split("\\s+");
				// str = sCurrentLine.split(" ");
				cities[j] = new City();
				// set coordinates to city
				cities[j].setX(Integer.parseInt(str[1]));
				cities[j].setY(Integer.parseInt(str[2]));
				// calculate weightMatrix
				for (int i = 0; i <= j; i++) {
					if (i == j) {
						weightMatrix[j][i] = 0;
					} else {
						weightMatrix[j][i] = weightMatrix[i][j] = Math
								.sqrt(Math.pow((cities[j].getX() - cities[i].getX()), 2)
										+ Math.pow((cities[j].getY() - cities[i].getY()), 2));
					}
				}
			}
			vertexInCluster = new int[numberOfCluster][];
			sCurrentLine = br.readLine();
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			root = Integer.parseInt(str[1]);
			// System.out.println(root);
			int max = 0, index = 0;
			for (int i = 0; i < numberOfCluster; i++) {
				int arrayCluster;
				sCurrentLine = br.readLine();
				str = sCurrentLine.split(" ");
				Cluster c = new Cluster();
				vertexInCluster[i] = new int[str.length - 2];
				if (max < str.length - 2) {
					max = str.length - 2;
					index = i;
				}
				for (int j = 1; j < str.length - 1; j++) {
					arrayCluster = Integer.parseInt(str[j]);
					c.getCluster().add(arrayCluster);
					vertexInCluster[i][j - 1] = arrayCluster;
				}
				try {
					clusters.add(c);
				} catch (Exception e) {
					System.out.println(e);
				}
			}
			// System.out.println(clusters.size());
			maxGroupIndex = index;
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		// return numberOfCity; // return the number of City
	}
}