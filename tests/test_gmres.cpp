// Small known systems exercise the same GMRES used by capacitance extraction.
// Faults injected at the operator boundary check true-residual verification and
// propagation of numerical, cancellation and out-of-core errors.
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <limits>

#include "SolveCapacitance.h"
#include "FasterCapGlobal.h"

static int failures = 0;

static void check(const char *name, bool passed)
{
    std::printf("%s %s\n", passed ? "ok    " : "FAILED", name);
    if(!passed) failures++;
}

class CSmallGmres : public CSolveCap
{
public:
    CSmallGmres(const CLin_Matrix &matrix, bool preconditioned = false)
        : m_clsMatrix(matrix), m_lProducts(0), m_lPreconditions(0),
          m_lFailProduct(0), m_lInvalidProduct(0), m_lPerturbProduct(0),
          m_bInvalidPreconditioner(false), m_bVariablePreconditioner(false)
    {
        m_clsGlobalVars.m_ucPrecondType = preconditioned ? AUTOREFINE_PRECOND_JACOBI : AUTOREFINE_PRECOND_NONE;
        m_bGmresConverged = true;
        check("allocate GMRES workspace", AllocateGmresMemory(matrix.num_rows()) == FC_NORMAL_END);
    }

    int SolveSystem(CLin_Vector &b, CLin_Vector &x, double tolerance, bool flexible = false, CLin_Vector *initial = NULL)
    {
        m_bGmresConverged = true;
        m_clsGlobalVars.m_dGmresTol = tolerance;
        if(initial != NULL) return gmresPrecondSFastAllX0(&b, &x, initial);
        return flexible ? gmresFlexPrecondSFastAll(&b, &x, tolerance) : gmresPrecondSFastAll(&b, &x, tolerance);
    }

    void Automatic(double tolerance) { m_clsGlobalVars.m_bAuto = true; m_clsGlobalVars.m_dAutoMaxErr = tolerance; }

    void CheckPhysicalCorrection(bool complex)
    {
        CAutoConductor dielectric, conductor;
        dielectric.m_bIsDiel = true;
        dielectric.m_ulLeafPanelNum = 1;
        conductor.m_bIsDiel = false;
        conductor.m_ulLeafPanelNum = 2;
        conductor.m_dSurfOutperm[0][0] = 10.0;
        conductor.m_dSurfOutperm[0][1] = complex ? -2.0 : 0.0;
        conductor.m_dSurfOutperm[1][0] = 1.0;
        conductor.m_dSurfOutperm[1][1] = 0.0;
        unsigned char materials[] = {0,0,1};
        m_clsMulthier.m_pucDielIndex = materials;
        m_clsMulthier.m_ulPanelNum[0] = 3;
        m_clsMulthier.m_stlConductors.push_back(&dielectric);
        m_clsMulthier.m_stlConductors.push_back(&conductor);
        m_clsGlobalVars.m_ucHasCmplxPerm = complex ? AUTOREFINE_CPLX_PERM : AUTOREFINE_REAL_PERM;
        CLin_Vector charge(complex ? 6 : 3, 0.0), correction(charge.size(), 0.0);
        charge[0] = 1e8; charge[1] = charge[2] = 1.0;
        correction[1] = complex ? 10.0/104.0 : 0.1;
        if(complex) correction[4] = 2.0/104.0;
        correction[2] = -1.0;
        double ratio = CSolveCap::GmresCorrectionRatio(correction, charge);
        double expected = 2.0/std::hypot(11.0, complex ? 2.0 : 0.0);
        check("weighted correction bounds cancellation despite large dielectric charge",
              std::fabs(ratio-expected) < 1e-14);
        // These synthetic material records and index bytes are stack-owned.
        m_clsMulthier.m_stlConductors.clear();
        m_clsMulthier.m_pucDielIndex = NULL;
        m_clsMulthier.m_ulPanelNum[0] = 0;
    }

    bool Converged() const { return m_bGmresConverged; }
    CLin_Matrix m_clsMatrix;
    long m_lProducts, m_lPreconditions, m_lFailProduct, m_lInvalidProduct, m_lPerturbProduct;
    bool m_bInvalidPreconditioner, m_bVariablePreconditioner;

protected:
    double GmresCorrectionRatio(const CLin_Vector &correction, const CLin_Vector &charge) override
    {
        // These synthetic systems have no conductor geometry; use their known
        // solution-vector norm as the observable. Physical tests cover the
        // production conductor/permittivity projection.
        double delta = 0.0, value = 0.0;
        for(unsigned long i=0; i<charge.size(); ++i) {
            delta = std::hypot(delta, correction[i]);
            value = std::hypot(value, charge[i]);
        }
        return delta == 0.0 ? 0.0 : delta/value;
    }

    int MultiplyGmresMatrix(CLin_Vector *result, CLin_Vector *x) override
    {
        m_lProducts++;
        if(m_lProducts == m_lFailProduct) return FC_CANNOT_GO_OOC;
        for(unsigned long i=0; i<x->size(); i++) {
            double sum = 0.0;
            for(unsigned long j=0; j<x->size(); j++) sum += m_clsMatrix[i][j] * (*x)[j];
            (*result)[i] = sum;
        }
        if(m_lProducts == m_lInvalidProduct) (*result)[0] = std::numeric_limits<double>::quiet_NaN();
        if(m_lProducts == m_lPerturbProduct) (*result)[0] += 0.25;
        return FC_NORMAL_END;
    }

    void PreconditionGmresVector(CLin_Vector *result, CLin_Vector *x) override
    {
        m_lPreconditions++;
        for(unsigned long i=0; i<x->size(); i++) {
            double scale = m_bVariablePreconditioner && m_lPreconditions % 2 ? 0.5 + 0.125 * i : 1.0;
            (*result)[i] = scale * (*x)[i] / m_clsMatrix[i][i];
        }
        if(m_bInvalidPreconditioner) (*result)[0] = std::numeric_limits<double>::infinity();
    }
};

static CLin_Matrix identity(unsigned long size)
{
    CLin_Matrix matrix(size, size);
    matrix = 0.0;
    for(unsigned long i=0; i<size; i++) matrix[i][i] = 1.0;
    return matrix;
}

static double relativeResidual(const CLin_Matrix &matrix, const CLin_Vector &b, const CLin_Vector &x)
{
    // Independent direct residual, using long double accumulation for the small test systems.
    long double normb = 0.0, normr = 0.0;
    for(unsigned long i=0; i<b.size(); i++) {
        long double residual = b[i];
        for(unsigned long j=0; j<x.size(); j++) residual -= (long double)matrix[i][j] * x[j];
        normr += residual * residual;
        normb += (long double)b[i] * b[i];
    }
    return std::sqrt((double)(normr / normb));
}

static void simpleSystems(bool flexible)
{
    CLin_Matrix matrix = identity(1);
    CSmallGmres solver(matrix);
    CLin_Vector b(1, 3.0), x(1, -9.0);
    std::feclearexcept(FE_ALL_EXCEPT);
    int ret = solver.SolveSystem(b, x, 1e-12, flexible);
    check("identity reaches happy breakdown", ret == FC_NORMAL_END && solver.Converged() && std::fabs(x[0]-3.0) < 1e-14);
    check("happy breakdown avoids zero division", !(std::fetestexcept(FE_INVALID | FE_DIVBYZERO)));
    check("happy breakdown verifies actual residual", solver.m_lProducts == 2);

    CSmallGmres zero(matrix, true);
    zero.m_bInvalidPreconditioner = true;
    b[0] = 0.0;
    ret = zero.SolveSystem(b, x, 1e-12, flexible);
    check("zero RHS gives zero solution", ret == FC_NORMAL_END && zero.Converged() && x[0] == 0.0);
    check("zero RHS skips preconditioner", zero.m_lPreconditions == 0);

    matrix[0][0] = 2.0;
    CSmallGmres exact(matrix, true);
    b[0] = 6.0;
    std::feclearexcept(FE_ALL_EXCEPT);
    ret = exact.SolveSystem(b, x, 1e-12, flexible);
    check("exact initial preconditioner solution", ret == FC_NORMAL_END && exact.Converged() && x[0] == 3.0);
    check("initial solution avoids zero division", !(std::fetestexcept(FE_INVALID | FE_DIVBYZERO)));
    check("initial solution uses one true residual", exact.m_lProducts == 1);

    matrix[0][0] = 0.0;
    CSmallGmres singular(matrix);
    ret = singular.SolveSystem(b, x, 1e-12, flexible);
    check("zero operator fails convergence with finite last solution", ret == FC_NORMAL_END && !singular.Converged() && x[0] == 0.0);
}

static void knownSystem(bool preconditioned, bool flexible)
{
    CLin_Matrix matrix(3, 3);
    const double entries[3][3] = {{4.0, 1.0, -1.0}, {2.0, 7.0, 1.0}, {1.0, -3.0, 12.0}};
    for(unsigned long i=0; i<3; i++) for(unsigned long j=0; j<3; j++) matrix[i][j] = entries[i][j];
    const double rhs[] = {1.5, -11.5, 13.0};
    const double expected[] = {1.0, -2.0, 0.5};
    CLin_Vector b(3, rhs), x(3);
    CSmallGmres solver(matrix, preconditioned);
    solver.m_bVariablePreconditioner = flexible;
    int ret = solver.SolveSystem(b, x, 1e-12, flexible);
    bool correct = ret == FC_NORMAL_END && solver.Converged();
    for(unsigned long i=0; i<3; i++) correct = correct && std::fabs(x[i]-expected[i]) < 1e-11;
    check("nonsymmetric system recovers known solution", correct);
    check("nonsymmetric system meets independently calculated residual", relativeResidual(matrix, b, x) < 1e-12);

    CSmallGmres initial(matrix, preconditioned);
    CLin_Vector x0(3, expected);
    std::feclearexcept(FE_ALL_EXCEPT);
    ret = initial.SolveSystem(b, x, 1e-12, false, &x0);
    check("explicit exact initial solution", ret == FC_NORMAL_END && initial.Converged() && initial.m_lProducts == 1);
    check("explicit initial solution avoids zero division", !(std::fetestexcept(FE_INVALID | FE_DIVBYZERO)));
}

static void continueAfterFailedResidual(bool flexible)
{
    CLin_Matrix matrix = identity(2);
    matrix[1][1] = 2.0;
    CLin_Vector b(2, -1.0), x(2);
    CSmallGmres solver(matrix);
    solver.m_lPerturbProduct = 2;
    // First recursive residual is about 0.316, but the perturbed true residual is
    // about 0.481. Rejecting it must preserve the Krylov state for the second step.
    int ret = solver.SolveSystem(b, x, 0.4, flexible);
    check("failed true residual continues iteration", ret == FC_NORMAL_END && solver.Converged()
          && solver.m_lProducts == 4 && std::fabs(x[0] + 1.0) < 1e-12 && std::fabs(x[1] + 0.5) < 1e-12);
    check("continued iteration meets independently calculated residual", relativeResidual(matrix, b, x) < 1e-12);
}

static void errors(bool flexible)
{
    CLin_Matrix matrix = identity(1);
    CLin_Vector b(1, 1.0), x(1);
    CSmallGmres inaccurate(matrix);
    inaccurate.m_lPerturbProduct = 2;
    int ret = inaccurate.SolveSystem(b, x, 1e-12, flexible);
    check("recursive residual cannot hide inaccurate final product", ret == FC_NORMAL_END && !inaccurate.Converged() && std::isfinite(x[0]));

    CSmallGmres failed(matrix);
    failed.m_lFailProduct = 2;
    ret = failed.SolveSystem(b, x, 1e-12, flexible);
    check("final residual propagates operator error", ret == FC_CANNOT_GO_OOC);

    CSmallGmres invalid(matrix);
    invalid.m_lInvalidProduct = 1;
    ret = invalid.SolveSystem(b, x, 1e-12, flexible);
    check("nonfinite operator rejected", ret == FC_INVALID_RESULT);

    CSmallGmres invalidFinal(matrix);
    invalidFinal.m_lInvalidProduct = 2;
    ret = invalidFinal.SolveSystem(b, x, 1e-12, flexible);
    check("nonfinite final residual rejected", ret == FC_INVALID_RESULT);

    CSmallGmres invalidPreconditioner(matrix, true);
    invalidPreconditioner.m_bInvalidPreconditioner = true;
    ret = invalidPreconditioner.SolveSystem(b, x, 1e-12, flexible);
    check("nonfinite preconditioner rejected", ret == FC_INVALID_RESULT && invalidPreconditioner.m_lProducts == 0);

    CSmallGmres invalidRhs(matrix);
    b[0] = std::numeric_limits<double>::quiet_NaN();
    ret = invalidRhs.SolveSystem(b, x, 1e-12, flexible);
    check("nonfinite RHS rejected", ret == FC_INVALID_RESULT);
    b[0] = 1.0;
    g_bFCContinue = false;
    ret = invalidRhs.SolveSystem(b, x, 1e-12, flexible);
    g_bFCContinue = true;
    check("cancellation preserved", ret == FC_USER_BREAK);
}

static void automaticAccuracy()
{
    // Tiny residual in the weak mode hides an order-one solution error.
    // Retrying against the original RHS at a merely smaller tolerance would
    // accept exactly the same wrong iterate. Solve the residual equation instead.
    CLin_Matrix matrix = identity(2);
    matrix[1][1] = 1e-12;
    for(bool flexible : {false, true}) {
        for(double scale : {1e-200, 1.0, 1e200}) {
            CSmallGmres automatic(matrix);
            automatic.Automatic(0.1);
            CLin_Vector b(2), x(2);
            b[0] = scale; b[1] = scale*1e-12;
            int ret = automatic.SolveSystem(b, x, 0.05, flexible);
            check("automatic solve resolves hidden weak mode", ret == FC_NORMAL_END && automatic.Converged()
                  && std::fabs(x[0]/scale-1.0) < 1e-12 && std::fabs(x[1]/scale-1.0) < 1e-12);
        }
        CLin_Vector b(2), x(2);
        b[0] = 1.0; b[1] = 1e-12;
        CSmallGmres manual(matrix);
        int ret = manual.SolveSystem(b, x, 0.05, flexible);
        check("manual tolerance keeps its existing residual contract", ret == FC_NORMAL_END
              && manual.Converged() && x[1] < 1e-10 && relativeResidual(matrix,b,x) < 0.05);

        CSmallGmres failed(matrix);
        failed.Automatic(0.1);
        failed.m_lFailProduct = 4; // First correction's matrix product.
        ret = failed.SolveSystem(b,x,0.05,flexible);
        check("automatic correction propagates operator errors", ret == FC_CANNOT_GO_OOC);
        CSmallGmres invalid(matrix);
        invalid.Automatic(0.1);
        invalid.m_lInvalidProduct = 4;
        ret = invalid.SolveSystem(b,x,0.05,flexible);
        check("automatic correction rejects nonfinite products", ret == FC_INVALID_RESULT);

        CLin_Matrix boundedMatrix = identity(2);
        boundedMatrix[1][1] = 2.0;
        CSmallGmres bounded(boundedMatrix);
        bounded.Automatic(1e-12);
        b = 1.0;
        // Deliberately demand far smaller forward error than each inner solve's
        // residual target: the refinement work limit must report failure.
        ret = bounded.SolveSystem(b,x,0.4,flexible);
        check("automatic correction work is bounded and cannot claim convergence",
              ret == FC_NORMAL_END && !bounded.Converged() && bounded.m_lProducts <= 60
              && std::isfinite(x[0]) && std::isfinite(x[1]));
    }
    CSmallGmres physical(identity(3));
    physical.CheckPhysicalCorrection(false);
    physical.CheckPhysicalCorrection(true);
}

int main()
{
    automaticAccuracy();
    simpleSystems(false);
    simpleSystems(true);
    knownSystem(false, false);
    knownSystem(true, false);
    knownSystem(true, true);
    continueAfterFailedResidual(false);
    continueAfterFailedResidual(true);
    errors(false);
    errors(true);

    const double magnitudes[] = {1e-200, 1e200};
    for(double magnitude : magnitudes) {
        CSmallGmres scaled(identity(1));
        CLin_Vector b(1, magnitude), x(1);
        int ret = scaled.SolveSystem(b, x, 1e-12);
        check("scaled RHS avoids norm underflow and overflow", ret == FC_NORMAL_END && scaled.Converged() && std::fabs(x[0]/magnitude-1.0) < 1e-12);
    }

    CLin_Matrix matrix = identity(2);
    matrix[1][1] = 0.0;
    CSmallGmres inconsistent(matrix);
    CLin_Vector b(2, 1.0), x(2);
    int ret = inconsistent.SolveSystem(b, x, 1e-12);
    check("inconsistent system retains finite last solution", ret == FC_NORMAL_END && !inconsistent.Converged() && std::isfinite(x[0]) && std::isfinite(x[1]));
    return failures ? 1 : 0;
}
