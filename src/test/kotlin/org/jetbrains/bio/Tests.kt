package org.jetbrains.bio

import kotlin.test.assertTrue

/**
 * In my sincerest opinion, [assertTrue] without a custom message should not be used
 * under any circumstances. There's nothing more annoying than seeing a test fail
 * with an oh-so-informative message of "Expected the value to be true." Yeah, and I
 * expected the message to help me identify the problem. Guess we're both disappointed.
 *
 * So, I've collected a few assert helper methods which cover popular use cases of [assertTrue].
 *
 * This object is partially copied from "bioinf-commons"; for the reasons see the comments to
 * https://github.com/JetBrains-Research/span/commit/fabb0b91827dd098dd3b96760b540b337018a6b2
 */
object Tests {

    fun assertIn(substring: String, fullString: String) {
        // Process Windows with different line separators correctly.
        substring.lines().forEach { s ->
            assertTrue(s in fullString, "Expected <$s> to be in <$fullString>.")
        }
    }

}