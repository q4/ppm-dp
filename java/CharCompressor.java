/* $Id$ */

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Reader;
import java.io.Writer;
import java.util.Set;
import java.util.Random;

/** Implements a generic compressor over streams of Unicode characters.
  * <p>This class interfaces a given data model (such as PPM) to an
  * arithmetic coder, forming a compressor over character sequences.
  * The data model should range over integers (-1..0xFFFF+)
  * representing Unicode code points, except for -1 which denotes
  * the additional EOF symbol.<br>
  * A static uniform base distribution over the permitted range of
  * integers is available as {@code CharCompressor.base}.</p>
  * <dl><dt><b>Note:</b></dt><dd>The range of characters is currently
  * restricted to the <i>Basic Multilingual Plane</i> (BMP), but this
  * should really be changed in future.  Currently, code points outside
  * the BMP range can be coded using surrogate pairs, but this is a
  * hack really and makes it harder for data models to do their job well.
  * </dd></dl>
  * @see ByteCompressor */
public class CharCompressor {

  long encoded = 0;
  long decoded = 0;

  /** Base distribution on unicode character points + EOF.
    * @see #eof */
  static final UniformInteger base = new UniformInteger(-1,0x0FFFF);

  /** EOF symbol. */
  final int eof = -1;

  /** Data model with adaptive coding over integers. */
  AdaptiveCode<Integer> model = null;

  /** Constructs a new character stream compressor with given model. */
  public CharCompressor(AdaptiveCode<Integer> model) {
    this.model = model;
  }

  /** Pretrains the model with a given training sequence. */
  public void train(Iterable<Character> src) {
    for (Character c : src) {
      model.learn((int) c);
    }
  }

  /** Compresses a sequence of integer code points.
    * @param src source of input code points
    * @param bw compressed bits are written here */
  public void compress(Iterable<Integer> src, BitWriter bw, int none)
                                                    throws IOException {
    Arith ac = new Arith();
    ac.start_encode(bw);
    //DebugEncoder debug = new DebugEncoder(ac);
    for (Integer cp : src) {
      model.encode(cp, ac);
      model.learn(cp);
      encoded++;
    }
    model.encode(eof, ac);
    ac.finish_encode();
    bw.close();
  }
  
  /** Compresses a sequence of characters.
    * @param src source of characters
    * @param bw compressed bits are written here
    * @param ac arithmetic coder to be used for compression */
  public void compress(Iterable<Character> src, BitWriter bw, Coder ac)
                                                    throws IOException {
    ac.start_encode(bw);
    //DebugEncoder debug = new DebugEncoder(ac);
    for (Character c : src) {
      model.encode((int) c, ac);
      model.learn((int) c);
      encoded++;
    }
    model.encode(eof, ac);
    ac.finish_encode();
    bw.close();
  }
  
  /** Decompresses a bitstream.
    * @param br source of compressed bits
    * @param w decoded characters are written here
    * @param ad arithmetic coder to be used for decompression */
  public void decompress(BitReader br, Writer w, Coder ad) throws IOException {
    //Decoder debug = new DebugDecoder(ad);
    ad.start_decode(br);
    Integer x = model.decode(ad);
    while (x != eof) {
      // FIXME: only works for BMP
      w.write((int) x);
      model.learn(x);
      decoded++;
      x = model.decode(ad);
    }
    ad.finish_decode();
    w.close();
  }
  
  /** Generates a character sequence from the data model.
    * This works by feeding random bits into the data model and
    * decompressing. The sequence can terminate early if the model
    * generates an EOF symbol.
    * @param rnd random number source
    * @param length maximum number of characters to generate
    * @param w generated characters are written here
    * @param ag arithmetic coder to be used for generation */
  public void generate(Random rnd, int length, Writer w, Coder ag) 
                                                    throws IOException {
    int sl = 0;
    RandomBitReader rbr = new RandomBitReader(rnd);
    ag.start_decode(rbr,true);  // start decoding, and pad with a 1-bit
    while (sl < length) {
      Integer x = model.decode(ag);
      if (x != eof) {
        w.write(x);
        model.learn(x);
        sl++;
      } else {
        break;
      }
    }
    ag.finish_decode();
    w.close();
  }

  /** Measures the information content (in bits) of <var>src</var> using
    * the data model's logMass method.
    * <b>Note</b>: The "logp" may differ from the actual compression
    * length.
    * @param src character sequence to be measured */
  public double measureLogMass(Iterable<Character> src) {
    double info = 0.0;
    for (Character c : src) {
      info += model.logMass((int) c);
      model.learn((int) c);
      encoded++;
    }
    // include EOF
    info += model.logMass(eof);
    return - info / Tools.LN2;
  }

  /** Measures the information content (in bits) of <var>src</var> using
    * the data model's {@code encode} method.
    * <b>Note</b>: this may change the state of the data model.
    * @param src character sequence to be measured
    * @param ac arithmetic coder to be used for compression measurement */
  public long measure(Iterable<Character> src, Coder ac) {
    BitCounter bc = new BitCounter();
    ac.start_encode(bc);
    for (Character c : src) {
      model.encode((int) c, ac);
      model.learn((int) c);
      encoded++;
    }
    model.encode(eof, ac);
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
  
  /** Pretrains the model with a given file.
    * @param tfnm training filename (or empty string, for <i>stdin</i>) */
  public void train(String tfnm) throws IOException {
    Iterable<Character> source = IOTools.charSequenceFromFile(tfnm);
    train(source);
  }

  /** Compresses specified input file to specified output file.
    * @param ifnm input filename (or empty string, for <i>stdin</i>)
    * @param ofnm output filename (or empty string, for <i>stdout</i>)
    * @param ac arithmetic coder to be used for compression */
  public void compress(String ifnm, String ofnm, Coder ac) throws IOException {
    Iterable<Character> source = IOTools.charSequenceFromFile(ifnm);
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
    Writer target = IOTools.getWriter(ofnm);
    decompress(source,target,ad);
  }
  
  /** Generate <var>n</var> chars to output file <var>ofnm</var>.
    * @param rnd random source
    * @param n number of symbols to generate
    * @param ofnm output filename (or empty string, for <i>stdout</i>)
    * @param ag arithmetic coder to be used for generation */
  public void generate(Random rnd, int n, String ofnm, Coder ag)
                                                   throws IOException {
    Writer target = IOTools.getWriter(ofnm);
    generate(rnd,n,target,ag);
    target.close();
  }
  
  /** Measures the information content (in bits) of <var>src</var> using
    * an arithmetic coder, and the data model's encode method.
    * Calls {@code measure(Iterable<Character>)}.
    * @param ifnm filename of file to be measured
    * @param ac arithmetic coder to be used for compression measurement
    * @see #measure(Iterable,Coder) */
  public long measure(String ifnm, Coder ac) throws IOException {
    Iterable<Character> source = IOTools.charSequenceFromFile(ifnm);
    long bits = measure(source,ac);
    return bits;
  }
  
  /** Measures the information content (in bits) of <var>src</var> using
    * the data model's logMass method.
    * <b>Note</b>: the "logp" may differ from the actual compression
    * length.
    * @param ifnm filename of file to be measured
    * @see #measureLogMass(Iterable) */
  public double measureLogMass(String ifnm) throws IOException {
    Iterable<Character> source = IOTools.charSequenceFromFile(ifnm);
    double bits = measureLogMass(source);
    return bits;
  }
  
  /** Returns a String identifier of this compressor instance.
    * Calls {@code model.toString()}. */
  public String toString() {
    return model.toString();
  }
  
  /** Debugs the compression model by running encoder + decoder
    * next to each other, comparing that they match up.
    * @param src source of input characters
    * @param modelcopy an exact copy of model, in the same state
    * @param verbose generate verbose output using
    *                DebugEncoder + DebugDecoder
    * @see RegionTracer
    * @see DebugEncoder
    * @see DebugDecoder */
  public void debug(Iterable<Character> src, AdaptiveCode<Integer> modelcopy,
                                         boolean verbose) throws IOException {
    RegionTracer rt = new RegionTracer();
    Encoder ec = rt;
    Decoder dc = rt;
    if (verbose) {
      ec = new DebugEncoder(rt);
      dc = new DebugDecoder(rt);
    }
    Integer x = null;
    for (Character c : src) {
      model.encode((int) c, ec);
      model.learn((int) c);
      x = modelcopy.decode(dc);
      modelcopy.learn(x);
    }
    model.encode(eof, ec);
    x = modelcopy.decode(dc);
    if (x == eof) {
      System.err.println("Found EOF where expected. Success.");
    } else {
      System.err.println("Found symbol "+x+" but expected EOF ("+eof+"). Failed.");
    }
  }
  
  /** Generate raw compression performance plot data,
    * writing to the supplied PrintStream.
    * The data comes in two tab-separated columns: the first
    * column is the offset (in characters) of the input data, the
    * second column the number of bits written so far.
    * @param src stream of source characters
    * @param ps PrintStream for plain text plot data
    * @param k thinning parameter. For k &gt; 0: samples every kth point.
    *          For k=0: samples with dynamically decreasing rate.
    * @param ac arithmetic coder to be used */
  public void mkplotdata(Iterable<Character> src, PrintStream ps, int k, Coder ac) {
    BitCounter bc = new BitCounter();
    ac.start_encode(bc);
    long pos = 0;  // position in input stream
    long prevpos = 0;  // position at previous symbol
    long prevbits = 0; // bitcount at previous symbol
    long lastbits = 0; // bitcount at last written symbol
    long lastpos = -1; // position at last written symbol
    long lastdelta = -1;
    int j = 0;
    for (Character c : src) {
      if (j == 0) {
        long delta = (bc.bitsWritten() - prevbits);
        // suppress points which lie on a straight line between two others
        if (delta != lastdelta) {
          if (pos > 0) {
            ps.println(prevpos+" \t"+prevbits);
            lastdelta = delta;
            lastbits = prevbits;
            lastpos = pos;
          } else {
            ps.println(pos+" \t"+bc.bitsWritten());
            lastdelta = delta;
            lastbits = bc.bitsWritten();
            lastpos = pos;
          }
          //ps.println(pos+" \t"+bc.bits());
        }
        prevpos  = pos;
        prevbits = bc.bitsWritten();
        if (k > 0) {
          j = k-1;              // static thinning
        } else {
          j = (int) (pos/1000); // dynamic thinning
        }
      } else {
        j--;
      }
      model.encode((int) c, ac);
      model.learn((int) c);
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


}
