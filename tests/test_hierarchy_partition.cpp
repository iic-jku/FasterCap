// Coincident centroids must still produce finite trees containing every input element.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <vector>

#include "Autorefine.h"

class CHierarchyProbe : public CAutoRefine
{
public:
    CAutoPanel *Build(std::vector<CAutoPanel *> &elements)
    {
        std::vector<CAutoPanel *> scratch(elements.size());
        C3DBBox bbox;
        for(CAutoPanel *element : elements) bbox += element->GetCentroid();
        m_pPanelArray1 = elements.data();
        m_pPanelArray2 = scratch.data();
        m_iLevel = -1;
        return RecurBuild3DSuperHier(0, elements.size(), &bbox);
    }
    CAutoSegment *Build(std::vector<CAutoSegment *> &elements)
    {
        std::vector<CAutoSegment *> scratch(elements.size());
        C3DBBox bbox;
        for(CAutoSegment *element : elements) bbox += element->GetCentroid();
        m_pSegmentArray1 = elements.data();
        m_pSegmentArray2 = scratch.data();
        m_iLevel = -1;
        return RecurBuild2DSuperHier(0, elements.size(), &bbox);
    }
};

static int failures = 0;

static void visit(CAutoElement *element, unsigned int level, unsigned int &depth,
                  std::set<CAutoElement *> &nodes, std::set<CAutoElement *> &leaves)
{
    if(element == NULL || !nodes.insert(element).second) {
        failures++;
        return;
    }
    depth = std::max(depth, level);
    if(element->IsLeaf()) leaves.insert(element);
    else {
        visit(element->m_pLeft, level + 1, depth, nodes, leaves);
        visit(element->m_pRight, level + 1, depth, nodes, leaves);
    }
}

template<class Element>
static void checkTree(const char *name, CHierarchyProbe &probe, std::vector<Element *> &elements)
{
    std::set<CAutoElement *> expected(elements.begin(), elements.end()), nodes, leaves;
    double measure = 0.0;
    for(Element *element : elements) measure += element->GetDimension();
    Element *root = probe.Build(elements);
    unsigned int depth = 0;
    visit(root, 0, depth, nodes, leaves);
    unsigned int bound = 0;
    for(size_t capacity = 1; capacity < elements.size(); capacity *= 2) bound++;
    bool passed = root != NULL && leaves == expected && nodes.size() == 2 * elements.size() - 1
        && depth <= bound + 1 && std::isfinite(root->GetDimension())
        && std::fabs(root->GetDimension() - measure) <= measure * 1e-12;
    std::printf("%s %s: %zu leaves, depth %u\n", passed ? "ok    " : "FAILED", name, leaves.size(), depth);
    if(!passed) failures++;
    if(root != NULL) root->DeleteElementsTree();
}

int main()
{
    CHierarchyProbe probe;
    const unsigned int counts[] = {2, 3, 17};
    for(unsigned int count : counts) {
        for(bool outlier : {false, true}) {
            std::vector<CAutoPanel *> triangles;
            std::vector<CAutoSegment *> segments;
            for(unsigned int i=0; i<count; i++) {
                const double shift = outlier && i + 1 == count ? 16.0 : 0.0;
                const double width = i + 1.0;
                CAutoPanel *triangle = new CAutoPanel;
                triangle->m_clsVertex[0].pos(shift - width, -1, 0);
                triangle->m_clsVertex[1].pos(shift + width, -1, 0);
                triangle->m_clsVertex[2].pos(shift, 2, 0);
                triangle->CalcPanelGeomPar();
                triangles.push_back(triangle);
                CAutoSegment *segment = new CAutoSegment;
                segment->m_clsVertex[0].pos(shift - width, 0);
                segment->m_clsVertex[1].pos(shift + width, 0);
                segment->CalcSegmentGeomPar();
                segments.push_back(segment);
            }
            checkTree(outlier ? "3D tied cluster" : "3D coincident centroids", probe, triangles);
            checkTree(outlier ? "2D tied cluster" : "2D coincident centroids", probe, segments);
        }
    }
    return failures ? 1 : 0;
}
