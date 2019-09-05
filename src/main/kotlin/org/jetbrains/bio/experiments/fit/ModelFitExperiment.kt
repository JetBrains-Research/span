package org.jetbrains.bio.experiments.fit

import org.apache.log4j.Logger
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiment.Experiment
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.Query
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.viktor.F64Array

/**
 * A generic experiment for evaluating classification models.
 *
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 12/09/13
 */
abstract class ModelFitExperiment<out Model : ClassificationModel, State : Any>(
        val genomeQuery: GenomeQuery,
        protected val dataQuery: Query<Chromosome, DataFrame>,
        protected val modelFitter: Fitter<Model>,
        private val modelClass: Class<out Model>,
        protected val availableStates: Array<State>)
    : Experiment("fit") {

    /** Secondary constructor to support effective genome query, see [SpanModelFitExperiment.createEffectiveQueries]. */
    internal constructor(
            queries: Pair<GenomeQuery, Query<Chromosome, DataFrame>>,
            modelFitter: Fitter<Model>,
            modelClass: Class<out Model>,
            availableStates: Array<State>
    ) : this(queries.first, queries.second, modelFitter, modelClass, availableStates)


    open val id: String
        get() {
            val model = modelClass.simpleName.replace(Regex("[^A-Z0-9]"), "").toLowerCase()
            require(availableStates.isNotEmpty())
            val states = availableStates[0].javaClass.simpleName.replace(Regex("[^A-Z0-9]"), "").toLowerCase()
            return reduceIds(listOf(genomeQuery.id, model, states, dataQuery.id))
        }


    fun getData(chromosome: Chromosome) = dataQuery.apply(chromosome)

    /**
     * Returns a data frame with the following columns:
     *
     *     offset  : int     0-based genomic offset
     *     state   : String  MLE state label
     *     ...     : double  state membership log-probability
     */
    abstract fun getStatesDataFrame(chromosome: Chromosome): DataFrame

    fun getStates(chromosome: Chromosome): List<State> {
        val statesMap = availableStates.associateBy { it.toString() }
        return getStatesDataFrame(chromosome).sliceAsObj<String>("state")
                .map { statesMap[it]!! }
    }

    fun getLogMemberships(chromosome: Chromosome): Map<State, F64Array> =
            getLogMemberships(getStatesDataFrame(chromosome))

    protected fun getLogMemberships(chromosomeStatesDF: DataFrame): Map<State, F64Array> =
            availableStates.associateBy({ it }) { chromosomeStatesDF.f64Array(it.toString()) }

    override fun doCalculations() {
        // IMPORTANT: since we have 2 types of experiments:
        // single model for all chromosomes and separate model for chromosome
        // getTrainedState() call can be blocking.
        // NO parallelism is allowed here!
        genomeQuery.get().forEach { chromosome ->
            getStatesDataFrame(chromosome)
        }
    }

    companion object {
        val LOG: Logger = Logger.getLogger(ModelFitExperiment::class.java)

    }
}
