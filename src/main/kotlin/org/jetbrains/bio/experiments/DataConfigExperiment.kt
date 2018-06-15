package org.jetbrains.bio.experiments

import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.jetbrains.bio.Logs
import org.jetbrains.bio.datasets.DataSet
import org.jetbrains.bio.util.*
import java.nio.file.Path

abstract class DataConfigExperiment(folder: String, val configuration: DataConfig) :
        Experiment("configs/${configuration.id}/$folder") {

    init {
        LOG.debug("Check configuration file")
        val configPath = experimentPath / "${configuration.id}.yaml"
        if (!configPath.isDirectory && configPath.isReadable) {
            val oldConfig: DataConfig?
            try {
                oldConfig = DataConfig.load(configPath)
            } catch (t: Throwable) {
                LOG.error("Failed to load config at $configPath", t)
                throw IllegalStateException("Failed to load config at $configPath", t)
            }
            check(configuration == oldConfig) {
                "Config file for $name $configPath changed.\nOLD config:\n$oldConfig\nNEW config:\n$configuration"
            }
            LOG.info("Config file already exists $name: $configPath")
        } else {
            saveConfig(configPath)
        }
    }

    private fun saveConfig(configPath: Path) {
        configPath.bufferedWriter().use { configuration.save(it) }
        LOG.info("Saved config file for $name: $configPath")
    }

    companion object {
        private val LOG = Logger.getLogger(DataConfigExperiment::class.java)

        fun loadDataConfig(input: String, quiet: Boolean=false): DataConfig {
            // Ensure that we have appender
            if (!quiet) {
                Logs.addConsoleAppender(Level.INFO)
            }
            // Let's try load from file
            val path = input.toPath()
            if (path.exists && path.isReadable) {
                if (!quiet) {
                    LOG.info("Loading data config from $path")
                }
                return DataConfig.load(path)
            }
            // Otherwise create by config
            val dataSet = DataSet.load(input)
            checkNotNull(dataSet) {
                "Cannot load data config for args: $input"
            }
            return dataSet!!.dataConfig
        }
    }
}