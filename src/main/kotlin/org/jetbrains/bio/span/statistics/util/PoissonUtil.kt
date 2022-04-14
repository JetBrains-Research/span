package org.jetbrains.bio.span.statistics.util

import org.jetbrains.bio.viktor.logAddExp
import kotlin.math.abs
import kotlin.math.ln

object PoissonUtil {
    /**
     * Poisson CDF evaluater for upper tail which allow calculation in log space.
     * @param k observation
     * @param lbd: Lambda
     * @return log(pvalue)
     *
     * See MACS2 sources Prob.pyx for original source code (log10PoissonCdfQLargeLambda):
     * ret = -lambda + \ln( \sum_{i=k+1}^{\inf} {lambda^i/i!} = -lambda + \ln( sum{ exp{ln(F)} } ), where F=lambda^m/m!
     * \ln{F(m)} = m*ln{lambda} - \sum_{x=1}^{m}\ln(x)
     * Calculate \ln( sum{exp{N} ) by logspace_add function
     */
    fun logPoissonCdf(k: Int, lbd: Double, maxM: Int = 10000, epsilon: Double = 1e-5): Double {
        require(lbd > 0) {
            "Lambda should be > 0, got $lbd"
        }
        var residue: Double
        var logX: Double
        val lnLbd = ln(lbd)
        // first residue
        val m = k + 1
        var sumLns = 0.0 // TODO[shpynov] this may be tabulated
        for (i in 1 until m + 1) {
            sumLns += ln(i.toDouble())
        }
        logX = m * lnLbd - sumLns
        residue = logX
        var logy: Double
        for (i in m + 1..maxM) { // Limit
            logy = logX + lnLbd - ln(i.toDouble())
            val preResidue = residue
            residue = preResidue logAddExp logy
            if (abs(preResidue - residue) < epsilon)
                break
            logX = logy
        }
        return residue - lbd
    }
}