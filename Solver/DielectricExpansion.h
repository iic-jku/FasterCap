// Coupled Cartesian expansion of the qualocation electric-field kernel.
// Source leaves represent point charges at their centroids; target leaves
// integrate the field over their actual triangles.
#ifndef FASTER_CAP_DIELECTRIC_EXPANSION_H
#define FASTER_CAP_DIELECTRIC_EXPANSION_H

#include "Geometry/Vector3D.h"
#include <cmath>

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

// Bounds use the physical descendants, not a superpanel's synthetic triangle.
// Qualocation source leaves are centroid points; target leaves are triangles.
class CDielectricExpansionBounds
{
public:
	void SetLeaf(const C3DVector_float vertices[3], C3DVector center)
	{
		m_dSourceRadius = m_dTargetRadius = 0.0;
		for(int i=0; i<3; ++i) {
			double distance = Mod(C3DVector(vertices[i])-center);
			if(distance > m_dTargetRadius) m_dTargetRadius = distance;
		}
	}

	void SetChildren(const CDielectricExpansionBounds &left, double leftShift,
		const CDielectricExpansionBounds &right, double rightShift)
	{
		m_dSourceRadius = std::fmax(left.m_dSourceRadius+leftShift, right.m_dSourceRadius+rightShift);
		m_dTargetRadius = std::fmax(left.m_dTargetRadius+leftShift, right.m_dTargetRadius+rightShift);
	}

	bool AcceptsSource(const CDielectricExpansionBounds &source, double distance) const
	{
		// Total degree two omits cubic terms. Keep both clusters inside one
		// quarter of the separation, well within the expansion's convergence ball.
		return m_dTargetRadius + source.m_dSourceRadius <= 0.25*distance;
	}

	double m_dSourceRadius, m_dTargetRadius;
};

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
		C3DVector r, double radial, bool expandTarget, bool sourceHasMoments = true)
	{
		const double x = r.x, y = r.y, z = r.z;
		const double invSquaredDistance = 1.0 / DotProd(r,r);
		double px = 0.0, py = 0.0, pz = 0.0, trace = 0.0;
		C3DVector secondTimesR(0,0,0);
		// Physical leaves are point sources. Their moments are identically zero;
		// skip those contractions without discarding neutral non-leaf groups.
		if(sourceHasMoments) {
			px = source.m_dDipole[0]; py = source.m_dDipole[1]; pz = source.m_dDipole[2];
			const double *s = source.m_dSecond;
			secondTimesR.x = s[0]*x + s[1]*y + s[2]*z;
			secondTimesR.y = s[1]*x + s[3]*y + s[4]*z;
			secondTimesR.z = s[2]*x + s[4]*y + s[5]*z;
			trace = s[0] + s[3] + s[5];
		}
		double dipoleDot = x*px + y*py + z*pz;
		double secondDot = x*secondTimesR.x + y*secondTimesR.y + z*secondTimesR.z;
		double alongR = (3*dipoleDot - 1.5*trace + 7.5*secondDot*invSquaredDistance) * invSquaredDistance;
		if(expandTarget) alongR += charge;
		m_clsField.x += radial*(x*alongR - px - 3*secondTimesR.x*invSquaredDistance);
		m_clsField.y += radial*(y*alongR - py - 3*secondTimesR.y*invSquaredDistance);
		m_clsField.z += radial*(z*alongR - pz - 3*secondTimesR.z*invSquaredDistance);
		if(!expandTarget) return;
		double identity = radial * (charge + 3*dipoleDot*invSquaredDistance);
		double outer = -radial * (3*charge + 15*dipoleDot*invSquaredDistance) * invSquaredDistance;
		double mixed = 3*radial*invSquaredDistance;
		// Store only independent Cartesian components: xx,xy,xz,yy,yz,zz.
		// Explicit components avoid indexed-vector branches in every linked pair.
		m_dGradient[0] += identity + outer*x*x + mixed*(x*px+x*px);
		m_dGradient[1] += outer*x*y + mixed*(x*py+y*px);
		m_dGradient[2] += outer*x*z + mixed*(x*pz+z*px);
		m_dGradient[3] += identity + outer*y*y + mixed*(y*py+y*py);
		m_dGradient[4] += outer*y*z + mixed*(y*pz+z*py);
		m_dGradient[5] += identity + outer*z*z + mixed*(z*pz+z*pz);
		double linear = -3*radial*charge*invSquaredDistance;
		double cubic = 15*radial*charge*invSquaredDistance*invSquaredDistance;
		// xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz.
		m_dHessian[0] += linear*(x+x+x) + cubic*x*x*x;
		m_dHessian[1] += linear*y + cubic*x*x*y;
		m_dHessian[2] += linear*z + cubic*x*x*z;
		m_dHessian[3] += linear*x + cubic*x*y*y;
		m_dHessian[4] += cubic*x*y*z;
		m_dHessian[5] += linear*x + cubic*x*z*z;
		m_dHessian[6] += linear*(y+y+y) + cubic*y*y*y;
		m_dHessian[7] += linear*z + cubic*y*y*z;
		m_dHessian[8] += linear*y + cubic*y*z*z;
		m_dHessian[9] += linear*(z+z+z) + cubic*z*z*z;
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
