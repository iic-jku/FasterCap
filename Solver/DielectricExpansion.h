// Coupled Cartesian expansion of the qualocation electric-field kernel.
// Source leaves represent point charges at their centroids; target leaves
// integrate the field over their actual triangles.
#ifndef FASTER_CAP_DIELECTRIC_EXPANSION_H
#define FASTER_CAP_DIELECTRIC_EXPANSION_H

#include "Geometry/Vector3D.h"

namespace DielectricExpansion {
inline unsigned int Index(int a, int b)
{
	static const unsigned int index[3][3] = {{0,1,2}, {1,3,4}, {2,4,5}};
	return index[a][b];
}

inline unsigned int Index(int a, int b, int c)
{
	static const unsigned int index[3][3][3] = {
		{{0,1,2}, {1,3,4}, {2,4,5}},
		{{1,3,4}, {3,6,7}, {4,7,8}},
		{{2,4,5}, {4,7,8}, {5,8,9}}
	};
	return index[a][b][c];
}
}

class CDielectricChargeMoments
{
public:
	void Clear()
	{
		for(int a=0; a<3; ++a) m_dDipole[a] = 0.0;
		for(int a=0; a<6; ++a) m_dSecond[a] = 0.0;
	}

	void AddChild(const CDielectricChargeMoments &child, double charge, C3DVector shift)
	{
		for(int a=0; a<3; ++a) {
			m_dDipole[a] += child.m_dDipole[a] + charge * shift[a];
			for(int b=a; b<3; ++b) {
				unsigned int ab = DielectricExpansion::Index(a,b);
				m_dSecond[ab] += child.m_dSecond[ab] + shift[a]*child.m_dDipole[b] +
					shift[b]*child.m_dDipole[a] + charge*shift[a]*shift[b];
			}
		}
	}

	double m_dDipole[3], m_dSecond[6];
};

class CDielectricFieldExpansion
{
public:
	void Clear()
	{
		m_clsField = C3DVector(0,0,0);
		for(int a=0; a<6; ++a) m_dGradient[a] = 0.0;
		for(int a=0; a<10; ++a) m_dHessian[a] = 0.0;
	}

	// Expand f(r+t-s) through total degree two in source and target shifts:
	// E = Q f - G p + H:S/2, G_local = Q G - H.p, H_local = Q H.
	// A leaf target already has the exact monopole integral in its scalar link;
	// it needs only the source-moment correction, without another local expansion.
	void AddSource(double charge, const CDielectricChargeMoments &source,
		C3DVector r, double radial, bool expandTarget)
	{
		double invSquaredDistance = 1.0 / DotProd(r,r);
		double gradient[3][3], hessian[3][3][3];
		for(int a=0; a<3; ++a) {
			if(expandTarget) m_clsField[a] += radial * charge * r[a];
			for(int b=0; b<3; ++b) {
				gradient[a][b] = radial * ((a == b ? 1.0 : 0.0) -
					3*r[a]*r[b]*invSquaredDistance);
				m_clsField[a] -= gradient[a][b] * source.m_dDipole[b];
				for(int c=0; c<3; ++c) {
					hessian[a][b][c] = radial * (-3 * ((a == b ? r[c] : 0.0) +
						(a == c ? r[b] : 0.0) + (b == c ? r[a] : 0.0))*invSquaredDistance +
						15*r[a]*r[b]*r[c]*invSquaredDistance*invSquaredDistance);
					m_clsField[a] += 0.5*hessian[a][b][c] *
						source.m_dSecond[DielectricExpansion::Index(b,c)];
				}
			}
		}
		if(!expandTarget) return;
		for(int a=0; a<3; ++a) {
			for(int b=a; b<3; ++b) {
				unsigned int ab = DielectricExpansion::Index(a,b);
				m_dGradient[ab] += charge * gradient[a][b];
				for(int c=0; c<3; ++c)
					m_dGradient[ab] -= hessian[a][b][c] * source.m_dDipole[c];
				for(int c=b; c<3; ++c)
					m_dHessian[DielectricExpansion::Index(a,b,c)] += charge * hessian[a][b][c];
			}
		}
	}

	void AddParent(const CDielectricFieldExpansion &parent, C3DVector shift)
	{
		m_clsField = m_clsField + parent.m_clsField;
		for(int a=0; a<3; ++a) {
			for(int b=0; b<3; ++b) {
				m_clsField[a] += parent.m_dGradient[DielectricExpansion::Index(a,b)] * shift[b];
				for(int c=0; c<3; ++c)
					m_clsField[a] += 0.5 * parent.m_dHessian[DielectricExpansion::Index(a,b,c)] * shift[b]*shift[c];
			}
			for(int b=a; b<3; ++b) {
				unsigned int ab = DielectricExpansion::Index(a,b);
				m_dGradient[ab] += parent.m_dGradient[ab];
				for(int c=0; c<3; ++c)
					m_dGradient[ab] += parent.m_dHessian[DielectricExpansion::Index(a,b,c)] * shift[c];
				for(int c=b; c<3; ++c) {
					unsigned int abc = DielectricExpansion::Index(a,b,c);
					m_dHessian[abc] += parent.m_dHessian[abc];
				}
			}
		}
	}

	double MeanNormalField(const C3DVector_float vertices[3], C3DVector center, C3DVector normal) const
	{
		C3DVector mean = m_clsField;
		// The triangle covariance about its centroid is sum(v_i v_i^T)/12.
		for(int vertex=0; vertex<3; ++vertex) {
			C3DVector shift = C3DVector(vertices[vertex]) - center;
			for(int a=0; a<3; ++a)
				for(int b=0; b<3; ++b)
					for(int c=0; c<3; ++c)
						mean[a] += m_dHessian[DielectricExpansion::Index(a,b,c)] * shift[b]*shift[c]/24.0;
		}
		return DotProd(mean,normal);
	}

	C3DVector m_clsField;
	double m_dGradient[6], m_dHessian[10];
};

#endif
