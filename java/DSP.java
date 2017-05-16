/* $Id$ */
import java.util.Random;

/** Interface common to discrete stochastic processes.
  * */
public interface DSP<X> {
  
  /** Advances the process and returns the next sample. */
  public X next(Random r);

}
