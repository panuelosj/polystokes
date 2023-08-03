#pragma once

#include "units.h"

// Utility for workin with corotational
//preconditioned conjugate gradient
template<typename PreconditionerSolver>
inline int pcg(Eigen::VectorXd& x, 
    const Eigen::SparseMatrix<SolveReal, Eigen::RowMajor>& A,
    const Eigen::VectorXd& b, 
    Eigen::VectorXd& r, 
    Eigen::VectorXd& z,
    Eigen::VectorXd& p, 
    Eigen::VectorXd& Ap, 
    PreconditionerSolver& pre,
    double& totalTimerAapply,
    double& totalTimerOther,
    double tol = 1e-1, 
    unsigned int num_itr = 100
)
{
    r = b - A * x;
    z = pre->solve(r);
    //z = r;
    p = z;
    double rsold = r.dot(z);
    double rsnew = 0.;
    double alpha = 0.;
    double beta = 0.;
    double rre = 0.;
    double xmag = 0.;

    // lets do something dumb and just run ONLY the application of A a whole bunch of times for no reason:
    std::chrono::time_point<std::chrono::high_resolution_clock> clockStart, clockEnd;
    double timerAapply = 0.;
    double timerOther = 0.;
    clockStart = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < 657; ++i) {
        Ap = A * p;
    }
    clockEnd = std::chrono::high_resolution_clock::now();
    timerAapply += std::chrono::duration<double, std::milli>(clockEnd - clockStart).count();
    totalTimerAapply = timerAapply;
    totalTimerOther = timerOther;

    for (unsigned int i = 0; i < num_itr; ++i) {
        Ap = A * p;

        alpha = rsold / (p.dot(Ap));
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r.dot(r);

        xmag = sqrt(x.dot(x));
        if (xmag > 0.) rre = rsnew / xmag;
        else rre = rsnew;

        //std::cout << i << ": " << rsnew << " / " << xmag << " = " << rre << std::endl;

        if (rre < tol) {
            return i;
        }
 
        //std::cout << "\t pre: z = " << z.dot(z) << ", r = " << r.dot(r) << std::endl;
        z = pre->solve(r);
        //std::cout << "\t post: z = " << z.dot(z) << ", r = " << r.dot(r) << std::endl;

        rsnew = r.dot(z);
        beta = rsnew / rsold;

        p = z + beta * p;
        rsold = rsnew;
    }

    return num_itr;
}


// Utility for workin with corotational
//preconditioned conjugate gradient
template<typename PreconditionerSolver>
inline int flex_pcg(Eigen::VectorXd& x, const Eigen::SparseMatrix<SolveReal, Eigen::RowMajor>& A,
    const Eigen::VectorXd& b, Eigen::VectorXd& r, Eigen::VectorXd& z,
    Eigen::VectorXd& p, Eigen::VectorXd& Ap, PreconditionerSolver& pre,
    double tol = 1e-1, unsigned int num_itr = 100)
{
    r = b - A * x;
    
    //std::cout << "\t pre: z = " << z.dot(z) << ", r = " << r.dot(r) << std::endl;
    z = pre->solve(r);
    //std::cout << "\t post: z = " << z.dot(z) << ", r = " << r.dot(r) << std::endl;

    Vector zold = z;
    Vector rold = r;
    p = z;
    double rsnew = 0.;
    double alpha = 0.;
    double beta = 0.;
    double rre = 0.;
    double xmag = 0.;

    for (unsigned int i = 0; i < num_itr; ++i) {
        Ap = A * p;
        alpha = rold.dot(zold) / (p.dot(Ap));
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r.dot(r);

        xmag = sqrt(x.dot(x));
        if (xmag > 0.) rre = rsnew / xmag;
        else rre = rsnew;

        //std::cout << i << ": " << rsnew << " / " << xmag << " = " << rre << std::endl;

        if (rre < tol) {
            return i;
        }
 
        //std::cout << "\t pre: z = " << z.dot(z) << ", r = " << r.dot(r) << std::endl;
        z = pre->solve(r);
        //std::cout << "\t post: z = " << z.dot(z) << ", r = " << r.dot(r) << std::endl;

        beta = r.dot(z - zold) / rold.dot(zold);
        zold = z;
        rold = r;

        //std::cout << "beta = " << beta << std::endl;
        p = z + beta * p;
        //std::cout << "p = " << p.dot(p) << std::endl;
    }
    return num_itr;
}

template<typename ApplyMatrix, typename PreconditionerSolver>
inline int bicgstab_external_matrix_A(Eigen::VectorXd& x,
    ApplyMatrix& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& r,
    Eigen::VectorXd& z,
    Eigen::VectorXd& p,
    Eigen::VectorXd& Ap,
    PreconditionerSolver& pre,
    double& totalTimerAapply,
    double& totalTimerOther,
    double& rre,
    double tol = 1e-1,
    unsigned int num_itr = 100
)
{
    r = b - A->apply(x);
    Eigen::VectorXd rhat = r;
    double rhoCurr = 1.;
    double rhoOld = 1.;
    double alpha = 1.;
    double beta = 0.;
    double omega = 1.;
    double xmag = 0.;
    double rsnew = 0.;
    p.resize(x.size());
    p.setZero();
    Eigen::VectorXd v(x.size());
    v.setZero();
    Eigen::VectorXd h(x.size());
    h.setZero();
    Eigen::VectorXd s(x.size());
    s.setZero();
    Eigen::VectorXd t(x.size());
    t.setZero();
    Eigen::VectorXd err(x.size());
    err.setZero();

    for (unsigned int i = 0; i < num_itr; ++i) {
        rhoOld = rhoCurr;
        rhoCurr = rhat.dot(r);
        beta = (rhoCurr / rhoOld) * (alpha / omega);
        p = r + beta * (p - omega * v);
        v = A->apply(p);
        alpha = rhoCurr / (rhat.dot(v));
        h = x + alpha * p;

        s = r - alpha * v;
        t = A->apply(s);
        omega = t.dot(s) / (t.dot(t));
        x = h + omega * s;

        // convergence check
        xmag = sqrt(x.dot(x));
        err = b - A->apply(x);
        rsnew = err.dot(err);
        rre = rsnew;
        if (sqrt(rsnew) / xmag < rre) rre = sqrt(rsnew) / xmag;
        if (rre < tol) {
            return i;
        }

        r = s - omega * t;
    }

    return num_itr;
}

template<typename ApplyMatrix, typename PreconditionerSolver>
inline int minres_external_matrix_A(Eigen::VectorXd& x,
    ApplyMatrix& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& r,
    Eigen::VectorXd& z,
    Eigen::VectorXd& p,
    Eigen::VectorXd& Ap,
    PreconditionerSolver& pre,
    double& totalTimerAapply,
    double& totalTimerOther,
    double& rre,
    double tol = 1e-1,
    unsigned int num_itr = 100
)
{
    Eigen::VectorXd p0, p1, p2;
    Eigen::VectorXd s0, s1, s2;
    double rsnew = 0.;
    double alpha = 0.;
    double beta1 = 0.;
    double beta2 = 0.;
    double xmag = 0.;

    r = b - A->apply(x);
    p0 = r;
    s0 = A->apply(p0);
    p1 = p0;
    s1 = s0;

    rre = 0.;

    for (unsigned int i = 0; i < num_itr; ++i) {
        p2 = p1; p1 = p0;
        s2 = s1; s1 = s0;

        alpha = r.dot(s1) / (s1.dot(s1));
        x += alpha * p1;
        r -= alpha * s1;

        xmag = sqrt(x.dot(x));
        rsnew = r.dot(r);
        rre = rsnew;
        if (sqrt(rsnew) / xmag < rre) rre = sqrt(rsnew) / xmag;
        if (rre < tol) {
            return i;
        }

        p0 = s1;
        s0 = A->apply(s1);
        beta1 = s0.dot(s1) / (s1.dot(s1));
        p0 -= beta1 * p1;
        s0 -= beta1 * s1;
        if (i > 1) {
            beta2 = s0.dot(s2) / (s2.dot(s2));
            p0 -= beta2 * p2;
            s0 -= beta2 * s2;
        }
    }

    return num_itr;
}


// Utility for workin with corotational
//preconditioned conjugate gradient
template<typename ApplyMatrix, typename PreconditionerSolver>
inline int pcg_external_matrix_A(Eigen::VectorXd& x,
    ApplyMatrix& A,
    const Eigen::VectorXd& b,
    Eigen::VectorXd& r,
    Eigen::VectorXd& z,
    Eigen::VectorXd& p,
    Eigen::VectorXd& Ap,
    PreconditionerSolver& pre,
    double& totalTimerAapply,
    double& totalTimerOther,
    double& rre,
    double tol = 1e-1,
    unsigned int num_itr = 100
)
{
    r = b - A->apply(x);
    z = pre->solve(r);
    //z = r;
    p = z;
    double rsold = r.dot(z);
    double rsnew = 0.;
    double alpha = 0.;
    double beta = 0.;
    double xmag = 0.;
    rre = 0.;

    // lets do something dumb and just run ONLY the application of A a whole bunch of times for no reason:
    /*
    std::chrono::time_point<std::chrono::high_resolution_clock> clockStart, clockEnd;
    double timerAapply = 0.;
    double timerOther = 0.;
    clockStart = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < 657; ++i) {
        Ap = A->apply(p);
    }
    clockEnd = std::chrono::high_resolution_clock::now();
    timerAapply += std::chrono::duration<double, std::milli>(clockEnd - clockStart).count();
    totalTimerAapply = timerAapply;
    totalTimerOther = timerOther;
    */

    for (unsigned int i = 0; i < num_itr; ++i) {
        Ap = A->apply(p);
        //Ap = A->apply(p) + (tol*0.1)*p; // diagonal hack

        alpha = rsold / (p.dot(Ap));
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r.dot(r);
        
        xmag = x.dot(x);
        rre = rsnew;
        if (rsnew / xmag < rre) rre = rsnew / xmag;
        if (rre < tol * tol) {
            rre = sqrt(rre);
            return i;
        }

        //std::cout << "\t pre: z = " << z.dot(z) << ", r = " << r.dot(r) << std::endl;
        z = pre->solve(r);
        //std::cout << "\t post: z = " << z.dot(z) << ", r = " << r.dot(r) << std::endl;

        rsnew = r.dot(z);
        beta = rsnew / rsold;

        p = z + beta * p;
        rsold = rsnew;
    }

    rre = sqrt(rre);
    return num_itr;
}
