/* JAVA Probability and Inference Tools
 * $Id$
 * Author: Christian Steinruecken */

import java.util.Random;
import java.util.Collection;

/** Sampler interface.
  * Methods implemented by all classes capable of generating samples.
  * */
public interface Sampler<X> {

  /** Returns a random sample from this distribution. */
  public abstract X sample(Random rnd);

  /** Samples <var>n</var> values and add them to collection <var>col</var>.
    * @see #sample(Random) */
  public void sample(Random rnd, int n, Collection<X> col);

}
