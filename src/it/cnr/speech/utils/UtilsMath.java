package it.cnr.speech.utils;

public class UtilsMath {

	public static int powerTwoApproximation(int n) {
		int y = 0;

		while (Math.pow(2, y) <= n) {
			y++;
		}
		y--; // lower approx
		return ((int) Math.pow(2, y));

	}
}
