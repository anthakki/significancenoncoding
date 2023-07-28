
//
// Multivariate normal estimator
//  for estimating dependency structure in the data
// 

public final class MVNEstimator {
	protected int _dim;
	protected double _sum0;
	protected double[] _sum1;
	protected double[] _sum2;

	protected int _sz(int n) {
		return n*(n+1)/2;
	}

	protected int _ix(int i, int j) {
		assert 0 <= i && i <= j;
		assert 0 <= j && j < _dim;
		return i + _sz(j);
	}

	public MVNEstimator(int dim) {
		_dim = dim;
		_sum0 = 0.;
		_sum1 = new double[_dim];
		_sum2 = new double[_sz(_dim)];
	}

	public void update(double[] x) {
		update(x, 1.);
	}

	public void update(double[] x, double w) {
		assert x.length == _dim;

		_sum0 += w;

		for (int i = 0; i < _dim; ++i)
			_sum1[i] += w * x[i];

		for (int j = 0; j < _dim; ++j)
			for (int i = 0; i <= j; ++i)
				_sum2[_ix(i, j)] += w * ( x[i] * x[j] );
	}

	public int dim() {
		return _dim;
	}

	public double size() {
		return _sum0;
	}

	public double mean(int i) {
		assert 0 <= i && i < _dim;
		return _sum1[i] / _sum0;
	}

	public double[] mean() {
		double[] m = new double[_dim];

		for (int i = 0; i < _dim; ++i)
			m[i] = mean(i);

		return m;
	}

	public double var(int i, double bessel) {
		return _cov(i, i, bessel);
	}

	public double var(int i) {
		return _cov(i, i, 1.);
	}

	protected double _cov(int i, int j, double bessel) {
		assert 0 <= i && i < j;
		assert 0 <= j && j < _dim;
		return ( _sum2[_ix(i, j)] - _sum1[i] / _sum0 * _sum1[j] ) / ( _sum0 - bessel );
	}

	public double cov(int i, int j, double bessel) {
		if (i <= j)
			return _cov(i, j, bessel);
		else
			return _cov(j, i, bessel);
	}

	public double cov(int i, int j) {
		return cov(i, j, 1.);
	}

	public double[][] cov(double bessel) {
		double[][] C = new double[_dim][_dim];
		double s = 1. / ( _sum0 - bessel );

		for (int j = 0; j < _dim; ++j) {
			double m_j = _sum1[j] / _sum0;

			for (int i = 0; i <= j; ++i) {
				double c = s * ( _sum2[_ix(i, j)] - _sum1[i] * m_j );
				C[i][j] = c;
				C[j][i] = c;
			}
		}

		return C;
	}

	public double[][] cov() {
		return cov(1.);
	}

	protected double _safeDiv(double x, double y) {
		// NB. compute x/y assuming |x| <= |y|
		if (y == 0.)
			return 0.;
		else
			return x / y;
	}

	public double cor(int i, int j, int diag) {
		if (i != j)
			return _safeDiv( cov(i, j, 0.), Math.sqrt( var(i, 0.) ) * Math.sqrt( var(j, 0.) ) );
		else if (diag == 0) {
			double s = Math.sqrt( var(i, 0.) );
			return _safeDiv( s * s, s * s );
		}
		else
			return 1.;
	}

	public double cor(int i, int j) {
		return cor(i, j, 1);
	}

	public double[][] cor(int diag) {
		double[][] C = cov(0.);

		for (int j = 0; j < _dim; ++j) {
			double C_jj = C[j][j] = Math.sqrt( C[j][j] );

			for (int i = 0; i < j; ++i) {
				double c = _safeDiv( C[j][i], C[i][i] * C_jj );
				C[j][i] = c;
				C[i][j] = c;
			}
		}

		if (diag == 0)
			for (int j = 0; j < _dim; ++j)
				C[j][j] = _safeDiv( C[j][j] * C[j][j], C[j][j] * C[j][j] );
		else
			for (int j = 0; j < _dim; ++j)
				C[j][j] = 1.;

		return C;
	}

	public double[][] cor() {
		return cor(1);
	}

	public int rank() {
		// NB. don't use cov(), it's rank-1 (rank) if mean is nonzero (zero)
		Jama.Matrix R = new Jama.Matrix(_dim, _dim);

		for (int j = 0; j < _dim; ++j)
			for (int i = 0; i <= j; ++i) {
				double r = _sum2[_ix(i, j)];
				R.set(i, j, r);
				R.set(j, i, r);
			}

		return R.svd().rank();
	}

	public double effSize() {
		// NB. singular values of an SPD matrix equal the eigenvalues
		double[] L = new Jama.Matrix(cor(0)).svd().getSingularValues();

		double s = 0., t = 0.;
		for (int i = 0; i < _dim; ++i) {
			s += Math.copySign( Math.sqrt(Math.abs( L[i] )), L[i] );
			t += L[i];
		}

		return s*s / t;
	}

} // MVNEstimator
