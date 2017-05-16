/* $Id$ */

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Random;
import java.util.Arrays;
import java.util.List;

/** Implements a generic bit stream compressor.
  * <p>This class interfaces a given data model (such as PPM) to an
  * arithmetic coder, forming a compressor over bit sequences.
  * The data model should range over the values Bit.ZERO and Bit.ONE.
  * The end of the bit stream is modelled using a separate model
  * over bit lengths.
  * A static base distribution over the bit values 0 and 1 
  * is available as {@code BitCompressor.base}.</p>
  * @see ByteCompressor */
public class BitCompressor {

  /** Low bit (zero). */
  static final Bit zero = Bit.ZERO;

  /** High bit (one). */
  static final Bit one  = Bit.ONE;

  /** EOF symbol. */
  static final Bit eof  = null;

  /** Base distribution on bit values (without EOF). */
  static final Bernoulli<Bit> base = new Bernoulli<Bit>(zero,one);
  
  /** Iterable collection of bit values + EOF. */
  static final List<Bit> values = Arrays.asList(new Bit[] { zero, one, eof });

  long encoded = 0;
  long decoded = 0;

  /** Data model with adaptive coding over bits. */
  AdaptiveCode<Bit> bitmodel = null;
  
  /** Model over sequence lengths. */
  Mass<Integer> eofmodel = null;

  /** Joint adaptive data model over bits + termination. */
  AdaptiveCode<Bit> model = null;

  /** Constructs a new Bit stream compressor with given model. */
  public BitCompressor(AdaptiveCode<Bit> bitmodel,
                       Mass<Integer> lengthmodel) {
    this.bitmodel = bitmodel;
    this.eofmodel = lengthmodel;
    this.model = new EOFModel<Bit>(bitmodel,eofmodel);
  }

  /** Pretrains the model with a given training sequence. */
  public void train(Iterable<Bit> src) {
    for (Bit b : src) { model.learn(b); }
  }

  /** Compresses a bit sequence.
    * @param src source of input bits
    * @param bw compressed bits are written here
    * @param ac arithmetic coder to be used for compression */
  public void compress(Iterable<Bit> src, BitWriter bw, Coder ac) throws IOException {
    ac.start_encode(bw);
    try {
      //DebugEncoder debug = new DebugEncoder(ac);
      for (Bit b : src) {
        model.encode(b, ac);
        model.learn(b);
        encoded++;
      }
      model.encode(eof, ac);
      ac.finish_encode();
      bw.close();
    }
    catch (ZeroMassException e) {
      throw new ZeroMassException("at bit "+encoded+" in input stream",e);
    }
  }
  
  /** Decompresses a bitstream.
    * @param br source of compressed bits
    * @param bw decoded bits are written here
    * @param ad arithmetic coder to be used for decompression */
  public void decompress(BitReader br, BitWriter bw, Coder ad) throws IOException {
    //Decoder debug = new DebugDecoder(ad);
    ad.start_decode(br);
    try {
      Bit x = model.decode(ad);
      while (x != eof) {
        bw.writeBit(x.byteValue());
        model.learn(x);
        decoded++;
        x = model.decode(ad);
      }
      ad.finish_decode();
      bw.close();
    }
    catch (ZeroMassException e) {
      throw new ZeroMassException("at bit "+decoded+" in output stream",e);
    }
  }
  
  /** Generates a bit sequence from the data model.
    * This works by feeding random bits into the data model and
    * decompressing. The sequence can terminate early if the model
    * generates an EOF symbol.
    * @param rnd random number source
    * @param length maximum number of bits to generate
    * @param bw generated bits are written here
    * @param ag arithmetic coder to be used for generation */
  public void generate(Random rnd, int length, BitWriter bw, Coder ag) 
                                                    throws IOException {
    int sl = 0;
    RandomBitReader rbr = new RandomBitReader(rnd);
    ag.start_decode(rbr,true);  // start decoding, and pad with a 1-bit
    while (sl < length) {
      Bit x = model.decode(ag);
      if (x != eof) {
        bw.writeBit(x);
        model.learn(x);
        sl++;
      } else {
        // IGNORE EOF
        //break;
      }
    }
    ag.finish_decode();
    bw.close();
  }
  
  /** Measures the information content (in bits) of <var>src</var> using
    * the data model's and length model's mass and logMass methods.
    * <b>Note</b>: this may change the state of the data model.
    * Also note that the "logp" may differ from the actual compression
    * length.
    * @param src bit sequence to be measured */
  public double measureLogMass(Iterable<Bit> src) {
    double info = 0.0;
    for (Bit b : src) {
      info += model.logMass(b);
      model.learn(b);
      encoded++;
    }
    info += model.logMass(eof);
    return - info / Tools.LN2;
  }
  
  /** Measures the information content (in bits) of <var>src</var> using
    * the data model's {@code encode} method.
    * <b>Note</b>: this changes the state of the data model.
    * @param src bit sequence to be measured
    * @param ac arithmetic coder to be used for compression measurement */
  public long measure(Iterable<Bit> src, Coder ac) {
    BitCounter bc = new BitCounter();
    ac.start_encode(bc);
    try {
      for (Bit x : src) {
        model.encode(x, ac);
        model.learn(x);
        encoded++;
      }
      model.encode(eof, ac);
    }
    catch (ZeroMassException e) {
      throw new ZeroMassException(e.toString()+" (at bit position "+encoded+")");
    }
    try {
      ac.finish_encode();
      bc.close();
    }
    catch (IOException e) {
      // should never happen...
      throw new RuntimeException(e);
    }
    return bc.bitsWritten();
  }

  /** Computes the current predictive distribution.
    * This is attempted first by calling the model's
    * {@code getPredictiveDistribution()} method.  If
    * the model does not implement this method (and
    * throws {@code UnsupportedOperationException} instead),
    * then a brute force copy is produced using the
    * model's {@code mass(Integer)} function.
    * @see AdaptiveCode#getPredictiveDistribution() */
  public Mass<Bit> getPredictiveDistribution() {
    Mass<Bit> predictive = null;
    try {
      predictive = model.getPredictiveDistribution();
    }
    catch (UnsupportedOperationException e) {
      /* This brute-force method will call model.mass(x) for all elements x. */
      Bit[] values = new Bit[] { zero, one, eof };
      double[] probs = new double[3];
      probs[0] = model.mass(zero);
      probs[1] = model.mass(one);
      probs[2] = model.mass(eof);
      predictive = new Discrete<Bit>(values,probs);
    }
    return predictive;
  }
  
  /** Prints the sequence of entropies (in bits) of each
    * conditional predictive distribution.
    * <b>Note</b>: this changes the state of the data model.
    * @see #getPredictiveDistribution()
    * @param src bit sequence
    * @param ps print stream */
  public void entropies(Iterable<Bit> src, PrintStream ps) {
    int pos = 0;
    Mass<Bit> dist;
    for (Bit b : src) {
      dist = getPredictiveDistribution();
      //ps.println(pos+" \t"+dist.mass(zero)+"\t"+dist.mass(one)+"\t"+dist.mass(eof));
      // compute and print entropy
      double entr = Tools.entropy(dist,values);
      ps.println(pos+" \t"+(entr / Tools.LN2));
      // update the model
      model.learn(b);
      pos++;
    }
    dist = getPredictiveDistribution();
    ps.println(pos+" \t"+(Tools.entropy(dist,values) / Tools.LN2));
    model.learn(eof);
    dist = getPredictiveDistribution();
    ps.println("# After EOF:");
    ps.println("--- \t"+(Tools.entropy(dist,values) / Tools.LN2));
  }
  
  /** Generates a deterministic adversarial (worst case)
    * bit sequence from the data model.
    * This works by picking the least predicted bit from the
    * predictive distribution every time.
    * When the two bits values are equiprobable, ZERO is chosen.
    * EOF is excluded.
    * @param length number of bits to generate
    * @param bw generated bits are written here */
  public void adverse(int length, BitWriter bw) throws IOException {
    int k = 0;
    while (k < length) {
      // get least predicted bit
      Bit lx = (model.mass(zero) > model.mass(one)) ? one : zero;
      bw.writeBit(lx);
      model.learn(lx);
      k++;
    }
    bw.close();
  }
  
  
  /** Generates a "friendly" deterministic best-case
    * bit sequence from the data model.
    * This works by picking a highest predicted bit from the
    * predictive distribution every time.  When the two bit
    * values are equiprobable, ZERO is chosen.
    * EOF is excluded.
    * @param length maximum number of bits to generate
    * @param bw generated bits are written here */
  public void friendly(int length, BitWriter bw) throws IOException {
    int k = 0;
    while (k < length) {
      // get least predicted bit
      Bit lx = (model.mass(zero) < model.mass(one)) ? one : zero;
      bw.writeBit(lx);
      model.learn(lx);
      k++;
    }
  }
  
  /** Pretrains the model with the content of a given file.
    * @param tfnm training filename (or empty string, for <i>stdin</i>) */
  public void train(String tfnm) throws IOException {
    Iterable<Bit> source = IOTools.bitSequenceFromFile(tfnm);
    train(source);
  }
  
  /** Compress input file <var>ifnm</var> to output file <var>ofnm</var>.
    * @param ifnm input filename (or empty string, for <i>stdin</i>)
    * @param ofnm output filename (or empty string, for <i>stdout</i>)
    * @param ac arithmetic coder to be used for compression */
  public void compress(String ifnm, String ofnm, Coder ac) throws IOException {
    Iterable<Bit> source = IOTools.bitSequenceFromFile(ifnm);
    BitWriter target = IOTools.getBitWriter(ofnm);
    compress(source,target,ac);
  }

  /** Decompress input file <var>ifnm</var> to output file <var>ofnm</var>.
    * @param ifnm input filename (or empty string, for <i>stdin</i>)
    * @param ofnm output filename (or empty string, for <i>stdout</i>)
    * @param ad arithmetic coder to be used for decompression */
  public void decompress(String ifnm, String ofnm, Coder ad)
                                                   throws IOException {
    BitReader source = IOTools.getBitReader(ifnm);
    BitWriter target = IOTools.getBitWriter(ofnm);
    decompress(source,target,ad);
  }
  
  /** Generate a sequence of bits to a specified output file.
    * @param rnd random source
    * @param n number of bits to generate
    * @param ofnm output filename (or empty string, for <i>stdout</i>)
    * @param ag arithmetic coder to be used for generation */
  public void generate(Random rnd, int n, String ofnm, Coder ag)
                                                       throws IOException {
    BitWriter target = IOTools.getBitWriter(ofnm);
    generate(rnd,n,target,ag);
    target.close();
  }
  
  /** Generate a sequence of deterministic adverse (worst case) bits
    * to a specified output file.
    * @param n number of bits to generate
    * @param ofnm output filename (or empty string, for <i>stdout</i>) */
  public void adverse(int n, String ofnm) throws IOException {
    BitWriter target = IOTools.getBitWriter(ofnm);
    adverse(n,target);
    target.close();
  }
  
  /** Generate a sequence of deterministic friendly (best case) bits
    * to a specified output file.
    * @param n number of bits to generate
    * @param ofnm output filename (or empty string, for <i>stdout</i>) */
  public void friendly(int n, String ofnm) throws IOException {
    BitWriter target = IOTools.getBitWriter(ofnm);
    friendly(n,target);
    target.close();
  }
  
  /** Measures the information content (in bits) of <var>src</var> using
    * the data model's {@code encode} method.
    * Calls {@code measure(Iterable<Bit>)}.
    * @param ifnm filename of file to be measured
    * @param ac arithmetic coder to be used for compression measurement
    * @see #measure(Iterable,Coder) */
  public long measure(String ifnm, Coder ac) throws IOException {
    Iterable<Bit> source = IOTools.bitSequenceFromFile(ifnm);
    return measure(source,ac);
  }
  
  /** Measures the information content (in bits) of <var>src</var> using
    * the data model's logMass method.
    * <b>Note</b>: this may change the state of the data model.
    * Also note that the "logp" may differ from the actual compression
    * length.
    * @param ifnm filename of file to be measured
    * @see #measureLogMass(Iterable) */
  public double measureLogMass(String ifnm) throws IOException {
    Iterable<Bit> source = IOTools.bitSequenceFromFile(ifnm);
    return measureLogMass(source);
  }

  /** Returns a String identifier of this compressor instance.
    * Calls {@code model.toString()}. */
  public String toString() {
    return model.toString();
  }
  
  /** Debugs the compression model by running encoder + decoder
    * next to each other, comparing that they match up.
    * @param src source of input bits
    * @param bitmodelcopy an exact copy of bit model
    * @param eofmodelcopy an exact copy of length model
    * @param verbose generate verbose output using
    *                DebugEncoder + DebugDecoder
    * @see RegionTracer
    * @see DebugEncoder
    * @see DebugDecoder */
  public void debug(Iterable<Bit> src,
                    AdaptiveCode<Bit> bitmodelcopy,
                    Mass<Integer> eofmodelcopy,
                    boolean verbose)  throws IOException {
    RegionTracer rt = new RegionTracer();
    Encoder ec = rt;
    Decoder dc = rt;
    if (verbose) {
      ec = new DebugEncoder(rt);
      dc = new DebugDecoder(rt);
    }
    EOFModel<Bit> modelcopy = new EOFModel<Bit>(bitmodelcopy,eofmodelcopy);
    Bit x;
    int pos=0;
    try {
      for (Bit b : src) {
        model.encode(b, ec);
        model.learn(b);
        x = modelcopy.decode(dc);
        modelcopy.learn(x);
        pos++;
      }
      pos = -1;
      model.encode(eof, ec);
      x = modelcopy.decode(dc);
      if (x == eof) {
        System.err.println("Found EOF where expected. Success.");
      } else {
        System.err.println("Found bit "+x+" but expected EOF. Failed.");
      }
    }
    catch (RuntimeException e) {
      System.err.println("Exception occurred at bit sequence position "+pos+":");
      throw e;
    }
  }

  /** Generate raw compression performance plot data,
    * writing to the supplied PrintStream.
    * The data comes in two tab-separated columns: the first
    * column is the offset (in bits) of the input data, the
    * second column the number of bits written so far.
    * @param src stream of source bits
    * @param ps PrintStream for plain text plot data
    * @param k thinning parameter. For k &gt; 0: samples every kth point.
    *          For k=0: samples with dynamically decreasing rate.
    * @param ac arithmetic coder to be used */
  public void mkplotdata(Iterable<Bit> src, PrintStream ps, int k, Coder ac) {
    BitCounter bc = new BitCounter();
    ac.start_encode(bc);
    long pos = 0;  // position in input stream
    long prevpos = 0;  // position at previous symbol
    long prevbits = 0; // bitcount at previous symbol
    long lastbits = 0; // bitcount at last written symbol
    long lastpos = -1; // position at last written symbol
    long lastdelta = -1;
    int j = 0;
    for (Bit b : src) {
      if (j == 0) {
        long delta = (bc.bitsWritten() - prevbits);
        ps.println(pos+" \t"+bc.bitsWritten());
        prevpos  = pos;
        prevbits = bc.bitsWritten();
        if (k > 0) {
          j = k-1;     // static thinning
        } else {
          if (pos < 200) { j = 0; } else
          if (pos < 500) { j = 4; } else
          if (pos < 1000) { j = 9; } else
          if (pos < 2000) { j = 19; } else
          if (pos < 5000) { j = 49; } else
          if (pos < 10000) { j = 99; } else
          if (pos < 20000) { j = 199; } else
          if (pos < 50000) { j = 499; } else
          if (pos < 100000) { j = 999; } else
          if (pos < 200000) { j = 1999; } else
          if (pos < 500000) { j = 4999; } else
          if (pos < 1000000) { j = 9999; } else
          if (pos < 2000000) { j = 19999; } else
          if (pos < 5000000) { j = 49999; } else
          if (pos < 10000000) { j = 99999; } else
          if (pos < 20000000) { j = 199999; } else
          if (pos < 50000000) { j = 499999; } else
          if (pos < 100000000) { j = 999999; } else
          {
            int l = (int) Math.floor(Math.log10(pos));
            int s = (int) Math.floor(Math.pow(10.0, l-1));
            j = s - 1;
          }
          // j = (int) (pos/1000); // dynamic thinning
        }
      } else {
        j--;
      }
      model.encode(b, ac);
      model.learn(b);
      encoded++;
      pos++;
    }
    ps.println(pos+" \t"+bc.bitsWritten());
    ps.println("# after encoding EOF:");
    model.encode(eof, ac);
    ps.println(pos+" \t"+bc.bitsWritten());
    ps.println("# after flushing the arithmetic coder to finish:");
    try {
      ac.finish_encode();
      ps.println(pos+" \t"+bc.bitsWritten());
      bc.close();
    }
    catch (IOException e) {
      // should never happen...
      throw new RuntimeException(e);
    }
  }
  
  /** Generate raw compression performance plot data, writing
    * to plot output file <var>pfnm</var>.
    * @see #mkplotdata(Iterable,PrintStream,int,Coder) */
  public void mkplotdata(Iterable<Bit> src, String pfnm, int k, Coder ac)
                                                        throws IOException {
    PrintStream plotdata = new PrintStream(IOTools.getOutputStream(pfnm));
    mkplotdata(src,plotdata,k,ac);
  }
  
  /** Generate raw compression performance plot data for a given
    * input file, writing to a given PrintStream.
    * @param ifnm input filename (or empty string, for <i>stdin</i>)
    * @param ps plot output PrintStream
    * @param ac arithmetic coder to be used
    * @see #mkplotdata(Iterable,PrintStream,int,Coder) */
  public void mkplotdata(String ifnm, PrintStream ps, Coder ac)
                                                        throws IOException {
    Iterable<Bit> src = IOTools.bitSequenceFromFile(ifnm);
    mkplotdata(src,ps,0,ac);
  }
  
  /** Generate raw compression performance plot data for a given
    * input file, writing to a given plot output file.
    * @param ifnm input filename (or empty string, for <i>stdin</i>)
    * @param pfnm plot output filename (or empty string, for <i>stdout</i>)
    * @param ac arithmetic coder to be used
    * @see #mkplotdata(Iterable,PrintStream,int,Coder) */
  public void mkplotdata(String ifnm, String pfnm, Coder ac) throws IOException {
    Iterable<Bit> src = IOTools.bitSequenceFromFile(ifnm);
    mkplotdata(src,pfnm,0,ac);
  }


  /** Generate symbol by symbol logp plot data,
    * writing to the supplied PrintStream.
    * The data comes in three tab-separated columns: the first
    * column is the offset (in bits) of the input data, the
    * second column the logp of the symbol at that offset,
    * the third column is the symbol itself (0 or 1).
    * @param src stream of source bits
    * @param ps PrintStream for plain text plot data */
  public void plotlogp(Iterable<Bit> src, PrintStream ps) {
    long pos = 0;  // position in input stream
    double logmass = Double.NEGATIVE_INFINITY;
    for (Bit x : src) {
      logmass = model.logMass(x) / Tools.LN2;
      ps.println(pos+" \t"+logmass+" \t"+x);
      model.learn(x);
      pos++;
    }
    ps.println("# after encoding EOF:");
    logmass = model.logMass(eof) / Tools.LN2;
    ps.println(pos+" \t"+logmass+" \tEOF");
  }
  
  /** Check that the predictive probability distributions
    * sum to one.
    * The data comes in two tab-separated columns: the first
    * column is the offset (in bits) of the input data, the
    * second column is the symbol itself, and the third column
    * is the sum of the model's current predictive probabilities.
    * This sum should always equal 1, or be very close to 1.
    * @param src stream of source bits
    * @param ps PrintStream for diagnostic data */
  public void sumcheck(Iterable<Bit> src, PrintStream ps) {
    long pos = 0;  // position in input stream
    for (Bit x : src) {
      double sum = model.mass(zero) + model.mass(one) + model.mass(eof);
      ps.println(pos+" \t"+x+" \t"+sum);
      model.learn(x);
      pos++;
    }
    ps.println("# after encoding EOF:");
    double sum = model.mass(zero) + model.mass(one) + model.mass(eof);
    ps.println(pos+" \tEOF \t"+sum);
  }
  
  /** Generate symbol by symbol debug information,
    * writing to the supplied PrintStream.
    * The data comes in two tab-separated columns: the first
    * column is the offset (in bits) of the input data, the
    * second column contains the model specific debug information.
    * @param src stream of source bits
    * @param ps PrintStream for diagnostic data */
  public void plotstate(Iterable<Bit> src, PrintStream ps) {
    long pos = 0;  // position in input stream
    java.lang.reflect.Method method = null;
    try {
      try {
        method = model.getClass().getMethod("getStateInfo");
      } catch (SecurityException e) {
      } catch (NoSuchMethodException e) {
      }
      for (Bit x : src) {
        ps.println(pos+" \t"+x+" \t"+method.invoke(model));
        model.learn(x);
        pos++;
      }
      ps.println("# before learning EOF:");
      ps.println(pos+" \tEOF \t"+method.invoke(model));
      model.learn(eof);
      ps.println("# after learning EOF:");
      ps.println(pos+" \t--- \t"+method.invoke(model));
    }
    catch (IllegalAccessException e) {
      throw new RuntimeException(e);
    }
    catch (java.lang.reflect.InvocationTargetException e) {
      throw new RuntimeException(e);
    }
  }
  
  /** Generate logp plot data, writing to output file <var>pfnm</var>.
    * @see #plotlogp(Iterable,PrintStream) */
  public void plotlogp(Iterable<Bit> src, String pfnm)
                                                       throws IOException {
    PrintStream plotdata = new PrintStream(IOTools.getOutputStream(pfnm));
    plotlogp(src,plotdata);
  }
  
  /** Generate logp plot data for input file <var>ifnm</var>,
    * writing to output file <var>pfnm</var>.
    * @see #plotlogp(Iterable,PrintStream) */
  public void plotlogp(String ifnm, String pfnm)
                                                       throws IOException {
    Iterable<Bit> src = IOTools.bitSequenceFromFile(ifnm);
    plotlogp(src,pfnm);
  }
  

}
