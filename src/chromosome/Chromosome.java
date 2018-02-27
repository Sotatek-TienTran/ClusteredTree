package chromosome;

public class Chromosome {

	private double[] factorialCost;
	double cost;
	private int[] factorialRank = new int[2];
	double scalarFitness;
	int skillFactor;

	private double[][] edgesMatrix; // Ma tran canh

	public void setFactorialRank(int rank, int index) {
		this.factorialRank[index] = rank;
	}

	public double getCost() {
		return cost;
	}

	public void setCost(double cost) {
		this.cost = cost;
	}

	public double[] getFactorialCost() {
		return factorialCost;
	}

	public void setFactorialCost(double[] factorialCost) {
		this.factorialCost = factorialCost;
	}

	public int[] getFactorialRank() {
		return factorialRank;
	}

	public void setFactorialRank(int[] factorialRank) {
		this.factorialRank = factorialRank;
	}

	public double getScalarFitness() {
		return scalarFitness;
	}

	public void setScalarFitness(double scalarFitness) {
		this.scalarFitness = scalarFitness;
	}

	public int getSkillFactor() {
		return skillFactor;
	}

	public void setSkillFactor(int skillFactor) {
		this.skillFactor = skillFactor;
	}

	public double[][] getEdgesMatrix() {
		return edgesMatrix;
	}

	public void setEdgesMatrix(double[][] edgesMatrix) {
		this.edgesMatrix = edgesMatrix;
	}

}
