package mfo_clustered;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.util.LinkedList;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class DrawLines extends JComponent {

	/**
		 * 
		 */
	private static final long serialVersionUID = 1L;
	static JFrame testFrame = new JFrame();

	private static class Line {
		final int x1;
		final int y1;
		final int x2;
		final int y2;
		String city1;
		String city2;
		final Color color;

		public Line(int x1, int y1, int x2, int y2, Color color, String city1, String city2) {
			this.x1 = x1;
			this.y1 = y1;
			this.x2 = x2;
			this.y2 = y2;
			this.color = color;
			this.city1 = city1;
			this.city2 = city2;
		}
	}

	private final LinkedList<Line> lines = new LinkedList<Line>();

	public void addLine(int x1, int x2, int x3, int x4, String city1, String city2) {
		addLine(x1, x2, x3, x4, Color.black, city1, city2);
	}

	public void addLine(int x1, int x2, int x3, int x4, Color color, String city1, String city2) {
		lines.add(new Line(x1, x2, x3, x4, color, city1, city2));
		repaint();
	}

	@Override
	protected void paintComponent(Graphics g) {
		super.paintComponent(g);
		for (Line line : lines) {
			g.setColor(line.color);
			//g.drawString
			g.drawLine(line.x1, line.y1, line.x2, line.y2);
			g.drawString(line.city1, line.x1+3, line.y1+3);
			g.drawString(line.city2, line.x2+3, line.y2+3);
			g.fillOval(line.x1-5, line.y1-5, 10 , 10);
		}
	}

	public  void draw(int x1, int y1, int x2, int y2, String city1, String city2) {

		testFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		final DrawLines comp = new DrawLines();
		comp.setPreferredSize(new Dimension(1000, 1000));
		testFrame.getContentPane().add(comp, BorderLayout.CENTER);
		JPanel buttonsPanel = new JPanel();

		testFrame.getContentPane().add(buttonsPanel, BorderLayout.SOUTH);

		comp.addLine(x1 * 10, y1 * 10, x2 * 10, y2 * 10, city1, city2);
		

		testFrame.pack();
		testFrame.setVisible(true);
	}
}
