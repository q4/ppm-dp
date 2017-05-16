/* $Id$ */

import java.util.Collection;

/** EOFmodel combines an adaptive model over symbols and a
  * model over sequence lengths.
  * Sequence termination is implicitly represented by
  * adjoining an EOF symbol. */
public class EOFModel<X> implements AdaptiveCode<X> {

  public final X EOF;

  /** An adaptive data model over symbols (without EOF). */
  public AdaptiveCode<X> datamodel;

  /** A probability distribution over sequence lengths. */
  public Mass<Integer> lengthmodel;

  /** Current file position. */
  int length = 0;

  /** Cumulative mass of termination at past file positions.
    * This mass is used to compute the conditional probability
    * of termination at the current position, given that
    * the sequence didn't terminate earlier. */
  double pastmass = 0.0;

  /** Prior probability of termination at the current position. */
  double lengthmass;

  /** Probability that the sequence continues at the current
    * position (as opposed to terminating here). */
  double pmore;

  /** Constructs a new EOFModel from a data model, length model, and
    * a designated EOF symbol that's unused by the data model.
    * In many cases, <code>null</code> can be used for eof.
    * @param datamodel an adaptive model over symbols
    * @param lengthmodel a distribution over sequence lengths
    * @param eof an unused symbol for termination, e.g. <code>null</code>. */
  public EOFModel(AdaptiveCode<X> datamodel,
                  Mass<Integer> lengthmodel, X eof) {
    this.datamodel = datamodel;
    this.lengthmodel = lengthmodel;
    this.EOF = eof;
    this.length = 0;
    this.lengthmass = lengthmodel.mass(length);
    this.pmore = 1.0 - lengthmass;
    this.pastmass = 0.0;
  }
  
  /** Constructs a new EOFModel from a data model and a length model.
    * EOF is internally represented as a <code>null</code> pointer.
    * @param datamodel an adaptive model over symbols
    * @param lengthmodel a distribution over sequence lengths */
  public EOFModel(AdaptiveCode<X> datamodel,
                  Mass<Integer> lengthmodel) {
    this(datamodel,lengthmodel,null);
  }

  /** Returns a String description of this model. */
  public String toString() {
    return "EOFModel with data="+datamodel.toString()
                 +" and length="+lengthmodel.toString();
  }
  
  /** Adapts the model by incorporating a single observation <var>x</var>. */
  public void learn(X x) {
    datamodel.learn(x);
    pastmass += lengthmass;
    length++;
    lengthmass = lengthmodel.mass(length);
    pmore = 1.0 - (lengthmass / (1.0 - pastmass));
  }

  /** Returns the predictive probability mass of the supplied symbol,
    * or the probability of EOF (if null). */
  public double mass(X x) {
    if (x != null) {
      /* This case is the most common one, and should be fast. */
      return datamodel.mass(x) * pmore;
    } else {
      return (1.0 - pmore); // probability of termination
    }
  }

  /** Returns the predictive log probability mass of the supplied symbol. */
  public double logMass(X x) {
    if (x != null) {
      return datamodel.logMass(x) + Math.log(pmore);
    } else {
      return Math.log(1.0 - pmore); // log probability of termination
    }
  }

  /** Returns the current predictive distribution. */
  public Mass<X> getPredictiveDistribution() {
    throw new UnsupportedOperationException();
  }

  /** Encodes a symbol (or sequence termination).
    * Termination is encoded when <var>x</var>
    * equals EOFModel.EOF.
    * Otherwise, continuation is encoded, followed
    * by the symbol value of <var>x</var>. */
  public void encode(X x, Encoder ec) {
    // first encode if the sequence continues
    double prob = pmore;
    if (prob == 0.0) {
      prob += Double.MIN_VALUE;
    }
    if (prob == 1.0) {
      prob -= Double.MIN_VALUE;
    }
    Bernoulli<Boolean> t = Bernoulli.booleans(prob);
    if (x == EOF) {
      // encode termination
      t.encode(false,ec);
    } else {
      // encode continuation
      t.encode(true,ec);
      // now encode the symbol
      datamodel.encode(x,ec);
    }
  }
 
  /** Decodes a symbol (or sequence termination).
    * @return either a symbol of type X, or EOFModel.EOF */
  public X decode(Decoder dc) {
    // first check if the sequence continues
    Bernoulli<Boolean> t = Bernoulli.booleans(pmore);
    boolean continues = t.decode(dc);
    if (continues) {
      // decode symbol
      return datamodel.decode(dc);
    } else {
      return EOF;
    }
  }
  
  public void encode(X x, Collection<X> omit, Encoder ec) {
    throw new UnsupportedOperationException();
  }
  
  public X decode(Collection<X> omit, Decoder dc) {
    throw new UnsupportedOperationException();
  }
  

}
