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
		double dipoleDot = 0.0, secondDot = 0.0;
		C3DVector secondTimesR(0,0,0);
		for(int a=0; a<3; ++a) {
			dipoleDot += r[a] * source.m_dDipole[a];
			for(int b=0; b<3; ++b)
				secondTimesR[a] += source.m_dSecond[DielectricExpansion::Index(a,b)] * r[b];
			secondDot += r[a] * secondTimesR[a];
		}
		double trace = source.m_dSecond[0] + source.m_dSecond[3] + source.m_dSecond[5];
		// Contract the symmetric moments directly instead of constructing all
		// 9 gradient and 27 Hessian entries for every directed link.
		double alongR = (3*dipoleDot - 1.5*trace + 7.5*secondDot*invSquaredDistance) * invSquaredDistance;
		if(expandTarget) alongR += charge;
		for(int a=0; a<3; ++a) {
			m_clsField[a] += radial * (r[a]*alongR - source.m_dDipole[a] -
				3*secondTimesR[a]*invSquaredDistance);
		}
		if(!expandTarget) return;
		double identity = radial * (charge + 3*dipoleDot*invSquaredDistance);
		double outer = -radial * (3*charge + 15*dipoleDot*invSquaredDistance) * invSquaredDistance;
		double mixed = 3*radial*invSquaredDistance;
		double hessianLinear = -3*radial*charge*invSquaredDistance;
		double hessianCubic = 15*radial*charge*invSquaredDistance*invSquaredDistance;
		for(int a=0; a<3; ++a) {
			for(int b=a; b<3; ++b) {
				unsigned int ab = DielectricExpansion::Index(a,b);
				m_dGradient[ab] += (a == b ? identity : 0.0) + outer*r[a]*r[b] +
					mixed*(r[a]*source.m_dDipole[b] + r[b]*source.m_dDipole[a]);
				for(int c=b; c<3; ++c) {
					m_dHessian[DielectricExpansion::Index(a,b,c)] += hessianLinear *
						((a == b ? r[c] : 0.0) + (a == c ? r[b] : 0.0) + (b == c ? r[a] : 0.0)) +
						hessianCubic*r[a]*r[b]*r[c];
				}
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
		C3DVector shift[3], sum(0,0,0);
		for(int vertex=0; vertex<3; ++vertex) {
			shift[vertex] = C3DVector(vertices[vertex]) - center;
			sum += shift[vertex];
		}
		// The stored float center need not equal the triangle's exact centroid.
		// Its first moment is sum/3; its second moment is
		// (sum*sum^T + sum_i(shift_i*shift_i^T))/12.
		double halfSecond[6];
		for(int b=0; b<3; ++b) {
			for(int c=b; c<3; ++c) {
				double moment = sum[b]*sum[c];
				for(int vertex=0; vertex<3; ++vertex)
					moment += shift[vertex][b]*shift[vertex][c];
				halfSecond[DielectricExpansion::Index(b,c)] = moment / 24.0;
			}
		}
		for(int a=0; a<3; ++a) {
			for(int b=0; b<3; ++b) {
				mean[a] += m_dGradient[DielectricExpansion::Index(a,b)] * sum[b]/3.0;
				for(int c=0; c<3; ++c)
					mean[a] += m_dHessian[DielectricExpansion::Index(a,b,c)] *
						halfSecond[DielectricExpansion::Index(b,c)];
			}
		}
		return DotProd(mean,normal);
	}

	C3DVector m_clsField;
	double m_dGradient[6], m_dHessian[10];
};

#endif
