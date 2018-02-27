package chromosome;

import java.util.Comparator;

public class ChromosomeCmp implements Comparator<Chromosome> {

	@Override
	public int compare(Chromosome ind1, Chromosome ind2) {
		// TODO Auto-generated method stub
		return (ind1.getScalarFitness() < ind2.getScalarFitness() ? -1
				: ind1.getScalarFitness() < ind2.getScalarFitness() ? 1 : 0);
	}

	public static Comparator<Chromosome> compareByFactorialCost = new Comparator<Chromosome>() {
		public int compare(Chromosome other, Chromosome one) {
			return Double.compare(other.cost, one.cost);
		}
	};
	public static Comparator<Chromosome> compareByScalarFitness = new Comparator<Chromosome>() {
		public int compare(Chromosome one, Chromosome other) {
			return Double.compare(other.getScalarFitness(), one.getScalarFitness());
		}
	};
	public static Comparator<Chromosome> compareByCost = new Comparator<Chromosome>() {
		public int compare(Chromosome other, Chromosome one) {
			return Double.compare(other.getCost(), one.getCost());
		}
	};
}