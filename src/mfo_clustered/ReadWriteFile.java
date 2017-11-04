package mfo_clustered;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class ReadWriteFile {
	public static int numberOfCity = clusterReadFiles("4eil76-2x2.clt");
	public static double weightMatrix[][];
	public static ArrayList<Cluster> clusters;
	public static int numberOfCluster;
	public static int root;
	public static int[][] vertexInCluster;
	public static City[] cities;

	public static int clusterReadFiles(String fileName) {
		clusters = new ArrayList<Cluster>();
		BufferedReader br = null;
		int num_city = 0;
		try {
			String sCurrentLine = null;
			br = new BufferedReader(new FileReader(fileName));
			for (int j = 0; j < 3; j++) {
				sCurrentLine = br.readLine();
			}
			String[] str = sCurrentLine.split(": ");
			num_city = Integer.parseInt(str[1]);

			weightMatrix = new double[num_city][num_city];
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			numberOfCluster = Integer.parseInt(str[1]);
			sCurrentLine = br.readLine();
			sCurrentLine = br.readLine();
			

			// read the detail of the City
			cities = new City[num_city];
			for (int j = 0; j < num_city; j++) {
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
			for (int i = 0; i < numberOfCluster; i++) {
				int arrayCluster;
				sCurrentLine = br.readLine();
				str = sCurrentLine.split(" ");
				Cluster c = new Cluster();
				vertexInCluster[i] = new int[str.length - 2];
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
		return num_city; // return number of City
	}
}
