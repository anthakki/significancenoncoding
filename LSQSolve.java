
//
// Generalized least-squares solvers
//  these generalize Jama.Matrix.solve() for non-square & singular matrices
//

public final class LSQSolve {
	public static Jama.SingularValueDecomposition _svd(Jama.Matrix A) {
		// NB. Jama chokes on special cases..

		if (!(A.getRowDimension() > 0 && A.getColumnDimension() > 0))
			return new Jama.SingularValueDecomposition( new Jama.Matrix( 1, 1, 0. ) ) {
				// NB. warn: S is 1-by-1 even though 0 rank!!!

				public Jama.Matrix getU() {
					return new Jama.Matrix( A.getRowDimension()    , getS().getRowDimension()    , 0. );
				}

				public Jama.Matrix getV() {
					return new Jama.Matrix( A.getColumnDimension() , getS().getColumnDimension() , 0. );
				}
			};

		if (!(A.getRowDimension() < A.getColumnDimension()))
			return A.svd();
		else
			return new Jama.SingularValueDecomposition( A.transpose() ) {
				public Jama.Matrix getU() {
					return super.getV().transpose();
				}

				public Jama.Matrix getV() {
					return super.getU().transpose();
				}
			};
	}

	private static Jama.Matrix _times(Jama.Matrix A, double d[]) {
		int m = A.getRowDimension();
		int n = A.getColumnDimension();
		assert d.length == n;

		for (int j = 0; j < n; ++j)
			for (int i = 0; i < m; ++i)
				A.set( i, j, A.get(i, j) * d[j] );

		return A;
	}

	private static double[] _solve(double d[], int m, int n) {
		double max_d = 0.;
		for (int j = 0; j < d.length; ++j)
			max_d = Math.max( max_d, d[j] );

		double th_d = Math.max(m, n) * Math.ulp(max_d);
		for (int j = 0; j < d.length; ++j)
			d[j] = d[j] > th_d ? 1. / d[j] : 0.;

		return d;
	}

	public static Jama.Matrix solve(Jama.Matrix A, Jama.Matrix B) {
		return solve(A).times(B);
	}

	public static Jama.Matrix solve(Jama.Matrix A) {
		Jama.SingularValueDecomposition svd = _svd(A);
		return _times( svd.getV(),
			_solve( svd.getSingularValues(), A.getRowDimension(), A.getColumnDimension() ) ).times(
				svd.getU().transpose() );
	}

	public static double[][] solve(double A[][], double B[][]) {
		return solve(new Jama.Matrix(A), new Jama.Matrix(B)).getArray();
	}

	public static double[][] solve(double A[][]) {
		return solve(new Jama.Matrix(A)).getArray();
	}

	public static double[] solve(double A[][], double b[]) {
		return solve(new Jama.Matrix(A), new Jama.Matrix(b, b.length)).getRowPackedCopy();
	}

} // LSQSolver
