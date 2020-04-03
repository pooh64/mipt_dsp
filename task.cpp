#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

typedef std::complex<float> cmplx;

void dft(cmplx *in, cmplx *out, size_t sz)
{
	for (size_t k = 0; k < sz; ++k) {
		out[k] = 0;
		float arg = -2 * M_PI * k / sz;
		cmplx a = std::exp(cmplx(0, 1) * arg);
		for (size_t n = 0; n < sz; ++n)
			out[k] += ((cmplx) std::pow(a, n)) * in[n];
	}
}

void dft_inv(cmplx *in, cmplx *out, size_t sz)
{
	for (size_t k = 0; k < sz; ++k) {
		out[k] = 0;
		float arg = 2 * M_PI * k / sz;
		cmplx a = std::exp(cmplx(0, 1) * arg);
		for (size_t n = 0; n < sz; ++n)
			out[k] += ((cmplx) std::pow(a, n)) * in[n];
		out[k] /= sz;
	}
}

void dacf(cmplx *in, cmplx *out, size_t sz)
{
	for (size_t k = 0; k < sz; ++k) {
		out[k] = 0;
		for (size_t n = k; n < sz; ++n)
			out[k] += in[n] * std::conj(in[n-k]);
	}
}

void sig_real(cmplx *in, float *out, size_t sz)
{
	for (size_t k = 0; k < sz; ++k)
		out[k] = std::real(in[k]);
}

void sig_imag(cmplx *in, float *out, size_t sz)
{
	for (size_t k = 0; k < sz; ++k)
		out[k] = std::imag(in[k]);
}

void sig_abs(cmplx *in, float *out, size_t sz)
{
	for (size_t k = 0; k < sz; ++k)
		out[k] = std::abs(in[k]);
}

void sig_sqrt(float *arr, size_t sz)
{
	for (size_t k = 0; k < sz; ++k)
		arr[k] = std::sqrt(arr[k]);
}

void rectangle_impulse(cmplx *out, size_t sz)
{
	size_t n1, n2;
	n1 = sz / 2 - sz / 8;
	n2 = sz / 2 + sz / 8;
	for (size_t n = 0; n < n1; ++n)
		out[n] = 0;
	for (size_t n = n1; n < n2; ++n)
		out[n] = 1;
	for (size_t n = n2; n < sz; ++n)
		out[n] = 0;
}

void triangle_sequence(cmplx *out, size_t sz, float len)
{
	for (size_t n = 0; n < sz; ++n) {
		long disc = long(float(n) / len);
		if (disc % 1 == 0) {
			float x = float(n) - disc * len;
			x = x / len - 0.5;
			out[n] = cmplx(1 - 2*std::abs(x), 0);
		} else {
			out[n] = cmplx(0, 0);
		}
	}
}

void out_to_file(char const *path, float *v,
	size_t from, size_t to)
{
	std::ofstream file;
	file.open(path);
	for (size_t i = from; i < to; ++i)
		file << v[i] << "\n";
	file.close();
}

int main()
{
	size_t mult = 8;		/* mult = 8 was used */
	size_t sz = 1024 * mult;	/* Number of points in signal buf */
	float len = 64;		 	/* Single triangle impulse lenght */

	auto s_signal = new cmplx[sz];
	auto s_buf1   = new cmplx[sz];
	auto s_buf2   = new cmplx[sz];
	auto s_result = new float[sz];

	/* Get simple dft(s) */
	triangle_sequence(s_signal, sz, len);		// s = triangle seq.
	dft(s_signal, s_buf1, sz);			// Apply dft
	sig_abs(s_buf1, s_result, sz);			// Aplly abs
	out_to_file("./out/abs_dft.dat", s_result, 0, 128 * mult);

	/* Check dft_inv(dft(s)) ?= s */
	dft_inv(s_buf1, s_buf2, sz);			// dft_inv(dft))
	sig_abs(s_buf2, s_result, sz);			// Apply abs
	out_to_file("./out/dftinv_dft.dat", s_result, 0, 128 * mult);

	/* Main task */
	dacf(s_signal, s_buf1, sz);			// Apply dacf
	dft(s_buf1, s_buf2, sz);			// Apply dft
	sig_real(s_buf2, s_result, sz);			// Apply real
	sig_sqrt(s_result, sz);				// Apply sqrt
	out_to_file("./out/sqrt_dft_dacf.dat", s_result, 0, 128 * mult);

	sig_real(s_buf1, s_result, sz);			// Get Re(dacf)
	out_to_file("./out/real_dacf.dat", s_result, 0, 128 * mult);
	sig_imag(s_buf1, s_result, sz);			// Get Im(dacf)
	out_to_file("./out/imag_dacf.dat", s_result, 0, 128 * mult);

	return 0;
}
