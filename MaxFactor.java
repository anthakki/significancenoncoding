
//
// New max factor models
//  this implements a new optimizer for the max factors and some additional models
//

import java.lang.Math;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;

public class MaxFactor {

	private static interface Density {
		// compute log-likelihood for the density
		public double log_f(double x, double m, double v);

		// compute log-likelihood density, gradient, Hessian wrt v 
		public double[] log_f_v(double[] f_v, double x, double m, double v);

		// compute tail probability (p-value)
		public double p_right(double x, double m, double v);
	};

	public static final Density DTY_GAMMA = new Density() {
		public double log_f(double x, double m, double v) {
			double b = m/v;
			double a = m*b;

			return a*Math.log(b) - Gamma.logGamma(a) + (a-1)*Math.log(x) - b*x;
		}

		public double[] log_f_v(double[] f_v, double x, double m, double v) {
			double b = m/v;
			double gb = -b/v;
			double hb = -2*gb/v;

			double a = m*b;
			// double ga = m*gb;
			// double ha = m*hb;

			double log_b = Math.log(b);
			double lgamma_a = Gamma.logGamma(a);
			double digamma_a = Gamma.digamma(a);
			double log_x = Math.log(x);

			double f = a*log_b - lgamma_a + (a-1)*log_x - b*x;
			double gf = gb*( m*( log_b + 1 - digamma_a + log_x ) - x );
			double hf = hb*( m*( log_b + 1.5 - digamma_a - .5*a*Gamma.trigamma(a) + log_x ) - x );

			f_v[0] = f;
			f_v[1] = gf;
			f_v[2] = hf;

			return f_v;
		}

		public double p_right(double x, double m, double v) {
			double b = m/v;
			double a = m*b;

			return Gamma.regularizedGammaQ(a, b*x);
		}

	}; // DTY_GAMMA

	public static final Density DTY_LOGT = new Density() {
		public static final double nu = 1.; // Cauchy

		public double log_f(double x, double m, double v) {
			double lm = Math.log( m*m / Math.sqrt( m*m + v ) );
			double ls = Math.sqrt( Math.log( 1. + v / (m*m) ) );

			double t = ( Math.log(x) - lm ) / ls;

			return -Math.log( ls*x ) +   // NB. dt/dx
				/* Gamma.logGamma( (nu+1)/2 ) - Gamma.logGamma( nu/2 ) - .5*Math.log( nu*Math.PI ) + */   // NB. constant
				-(nu+1)/2 * Math.log( 1. + t*t / nu );
		}

		public double[] log_f_v(double[] f_v, double x, double m, double v) {
			double lm = Math.log( m*m / Math.sqrt( m*m + v ) );
			double glm = -.5 / ( m*m + v );

			double ls = Math.sqrt( Math.log( 1 + v / ( m*m ) ) );
			double gls = -glm / ls;

			double t = ( Math.log(x) - lm ) / ls;
			double q = 1 + t*t / nu;

			double f = -Math.log( ls*x ) +
				/* Gamma.logGamma( (nu+1)/2 ) - Gamma.logGamma( nu/2 ) - .5*Math.log( nu*Math.PI ) + */
				-(nu+1)/2 * Math.log(q);

			// NB. ( t/ls )*( t/ls-1 ) = ( log(x)-lm )*( log(x)-(lm+ls*ls) )/( ls*ls )

			// TODO: this is correct, but I don't grok all of it.. would need some work..

			double gf = -gls/ls * ( 1 - (nu+1)/nu * ( t*ls )*( t/ls-1 )/q );
			double hf = -2 * gls/ls * gls * (( 1/ls + ls * ( 1 +
				- (nu+1)/nu * ( ( -.5 + ( t*ls )*( t/ls-1 ) + ( t/ls )*( t/ls ) )*q + ( t/ls-1 )*( t/ls-1 ) ) / ( q*q )
				)) );

			f_v[0] = f;
			f_v[1] = gf;
			f_v[2] = hf;

			return f_v;
		}

		public double p_right(double x, double m, double v) {
			double lm = Math.log( m*m / Math.sqrt( m*m + v ) );
			double ls = Math.sqrt( Math.log( 1 + v / ( m*m ) ) );

			double t = ( Math.log(x) - lm ) / ls;

			// NB. F(T<t) = 1-1/2*I_{nu/(nu+t*t)}(nu/2,1/2) for t>0

			if (!(t < 0.))
				return      .5 * Beta.regularizedBeta( nu / ( nu + t*t ), .5*nu, .5 );
			else
				return 1. - .5 * Beta.regularizedBeta( nu / ( nu + t*t ), .5*nu, .5 );
		}

	}; // DTY_LOGT

	public static interface VarianceModel {
		public double[] init_F();

		public double v(double m, double[] F);

		public double v_F(double[] v_F, double m, double[] F);
	};

	public static final VarianceModel VAR_P12 = new VarianceModel() {
		public double[] init_F() {
			return new double[]{ 0.55, 0.51 };
		}

		public double v(double m, double[] F) {
			return m * ( F[0] + m * ( F[1] ) );
		}

		public double v_F(double[] v_F, double m, double[] F) {
			v_F[0] = m;
			v_F[1] = m*m;

			return v(m, F);
		}

	}; // VAR_P12

	public static final VarianceModel VAR_P012 = new VarianceModel() {
		public double[] init_F() {
			double[] F = VAR_P12.init_F();
			return new double[]{ 1e-7, F[0], F[1] };
		}

		public double v(double m, double[] F) {
			return F[0] + m * ( F[1] + m * ( F[2] ) );
		}

		public double v_F(double[] v_F, double m, double[] F) {
			v_F[0] = 1.;
			v_F[1] = m;
			v_F[2] = m*m;

			return v(m, F);
		}

	}; // VAR_P012

	private Density _density;
	private VarianceModel _varModel;

	MaxFactor(Density density, VarianceModel varModel) {
		_density = density;
		_varModel = varModel;
	}

	public double[] solve(double[] x, double[] m, int verbose) {
		final int max_iter = 100;
		final double obj_tol = 1e-7;

		double[] F = _varModel.init_F();
		double[] new_F = new double[ F.length ];

		double[] dF = new double[ F.length ];
		double[] g = new double[ F.length ];
		double[] H = new double[ (F.length+1)*F.length/2 ];

		if (verbose >= 3)
			System.err.printf("%s.solve(density = %s, varModel = %s)\n",
				this.getClass().getName(),
				_density == DTY_GAMMA ? "DTY_GAMMA" : _density == DTY_LOGT ? "DTY_LOGT" : "(unknown)",
				_varModel == VAR_P12 ? "VAR_P12" : _varModel == VAR_P012 ? "VAR_P012" : "(unknown)");

		for (int iter = 0; iter < max_iter; ++iter) {
			double obj = _computeObjGrad(g, H, x, m, F);

			if (verbose >= 2)
				System.err.printf("f=%g g=%s H=%s\n", obj, java.util.Arrays.toString(g), java.util.Arrays.toString(H));

			_solveLog(dF, H, g, F);
			double a = 1.;

			double new_obj = obj;
			while (a > 0. && !( (new_obj = _computeObj(x, m, _mulLogDelta(new_F, F, a, dF))) >= obj ))
				a = .5 *  a;

			if (verbose >= 1)
				System.err.printf("iter=%d F=%s obj=%g (%+g)\n", iter, java.util.Arrays.toString(new_F), new_obj, new_obj - obj);

			if (!(new_obj - obj > obj_tol)) {
				if (verbose >= 1)
					System.err.printf("early exit\n");
				break;
			}

			{ double[] temp = F;
			F = new_F;
			new_F = temp; }
		}

		return F;
	}

	public double[] solve(double[] x, double[] m) {
		return solve(x, m, 0);
	}

	public double eval(double x, double m, double[] F) {
		return _computeTailP(x, m, F);
	}

	private double _computeObj(double[] x, double[] m, double[] F) {
		double logLike = 0.;

		for (int i = 0; i < x.length; ++i)
			logLike += _density.log_f( x[i], m[i], _varModel.v( m[i], F ) );

		return logLike;
	}

	private double _computeObjGrad(double[] g, double[] H, double[] x, double[] m, double[] F) {
		double[] f_v = { Double.NaN, Double.NaN, Double.NaN };   // ( d^0/dv^0 , d^1/dv^1 , d^2/dv^2 ) f
		double[] v_F = new double[ F.length ];   // d/dF v

		double logLike = 0.;

		for (int k = 0; k < g.length; ++k)
			g[k] = 0.;
		for (int k = 0; k < H.length; ++k)
			H[k] = 0.;

		for (int i = 0; i < x.length; ++i) {
			double v = _varModel.v_F( v_F, m[i], F );
			_density.log_f_v( f_v, x[i], m[i], v );

			logLike += f_v[0];

			// d/dF f = ( d/dF v ) * ( d/dv f )
			// d/dF d/dF' f = ( d/dF v ) * ( d/dF' v ) * ( d/dv d/dv f ) + ( d/dF d/dF' v ) * ( d/dv f )

			for (int k = 0; k < g.length; ++k)
				g[k] += v_F[k] * f_v[1];
			for (int k = 0, k2 = 0; k2 < g.length; ++k2)
			for (int        k1 = 0; k1 < k2 + 1  ; ++k1, ++k)
				H[k] += v_F[k1] * v_F[k2] * f_v[2];
		}

		return logLike;
	}

	private static double[] _solve(double[] p, double[] H, double[] g) {
		Jama.Matrix H1 = new Jama.Matrix(g.length, g.length);

		for (int k = 0, k2 = 0; k2 < g.length; ++k2)
		for (int        k1 = 0; k1 < k2 + 1  ; ++k1, ++k) {
			H1.set(k1, k2, H[k]);
			H1.set(k2, k1, H[k]);
		}

		Jama.Matrix p1 = LSQSolve.solve(H1).times(new Jama.Matrix(g, g.length));

		double dot = 0.;
		for (int k = 0; k < p.length; ++k)
			dot += g[k] * p1.get(k, 0);

		double sign = dot > 0. ? +1. : -1.;   // positive gradient dir
		for (int k = 0; k < p.length; ++k)
			p[k] = sign * p1.get(k, 0);

		return p;
	}

	private static double[] _solveLog(double[] p, double[] H, double[] g, double[] x) {
		double[] g1 = new double[g.length];

		for (int k = 0; k < g.length; ++k)
			g1[k] = x[k] * g[k];

		Jama.Matrix H1 = new Jama.Matrix(g.length, g.length);

		for (int k = 0, k2 = 0; k2 < g.length; ++k2)
		for (int        k1 = 0; k1 < k2 + 1  ; ++k1, ++k) {
			double v = x[k1] * x[k2] * H[k];
			H1.set(k1, k2, v);
			H1.set(k2, k1, v);
		}

		for (int k = 0; k < g.length; ++k)
			H1.set(k, k, H1.get(k, k) + g1[k]);

		Jama.Matrix p1 = LSQSolve.solve(H1).times(new Jama.Matrix(g, g.length));

		double dot = 0.;
		for (int k = 0; k < p.length; ++k)
			dot += g[k] * p1.get(k, 0);

		double sign = dot > 0. ? +1. : -1.;   // positive gradient dir
		for (int k = 0; k < p.length; ++k)
			p[k] = sign * p1.get(k, 0);

		return p;
	}

	private static double _boundToPositive(double[] F, double[] dF) {
		double a = 1.;

		for (int i = 0; i < F.length; ++i)
			if (dF[i] < 0.) {
				double a1 = -F[i] / dF[i];
				if (a1 < a)
					a = a1;
			}

		return a;
	}

	private static double[] _addDelta(double[] new_x, double[] x, double a, double[] dx) {
		for (int i = 0; i < x.length; ++i)
			new_x[i] = x[i] + a * dx[i];

		return new_x;
	}

	private static double[] _mulLogDelta(double[] new_x, double[] x, double a, double[] dx) {
		for (int i = 0; i < x.length; ++i)
			new_x[i] = x[i] * Math.exp( a * dx[i] );

		return new_x;
	}

	private double _computeTailP(double x, double m, double v) {
		return _density.p_right(x, m, v);
	}

	private double _computeTailP(double x, double m, double[] F) {
		return _computeTailP(x, m, _varModel.v(m, F));
	}

} // MaxFactor
