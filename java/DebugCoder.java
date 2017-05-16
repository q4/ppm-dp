/* $Id$ */

import java.io.IOException;

/** A debug wrapper for arithmetic coders.
  * The wrapper methods in this class print debugging output
  * with TTY colour control codes to stderr. */
public class DebugCoder implements Coder {

  private Coder coder = null;
  long lastrange = -1;

  public DebugCoder(Coder c) {
    this.coder = c;
  }

  public String toString() {
    return coder.toString() + " + debugging output";
  }

  private void warn(String warning) {
    System.err.println("\033[41;37;1mAC.WARNING! "+warning+"\033[m");
  }

  public void start_encode(BitWriter bw) {
    System.err.println("\033[33mAC.start_encode("+bw+")\033[m");
    coder.start_encode(bw);
  }
  public void finish_encode() throws IOException {
    System.err.println("\033[33mAC.finish_encode()\033[m");
    coder.finish_encode();
  }

  public void start_decode(BitReader br) throws IOException {
    System.err.println("\033[33mAC.start_decode("+br+")\033[m");
    coder.start_decode(br);
  }
  public void start_decode(BitReader br, boolean pad) throws IOException {
    System.err.println("\033[33mAC.start_decode("+br+",pad="+pad+")\033[m");
    coder.start_decode(br,pad);
  }
  public void finish_decode() {
    System.err.println("\033[33mAC.finish_decode()\033[m");
    coder.finish_decode();
  }


  public void storeRegion(long l, long h, long t) {
    double logp = Math.log((double) (h-l) / t) / Tools.LN2;
    System.err.println("\033[35mAC.stor "+l+", "+h+", "+t
                      +" \t\t(+"+(-logp)+" bits)\033[m");
    coder.storeRegion(l,h,t);
    lastrange = -1; // discard last known range
  }

  public void storeRegion(long l, long h) {
    if (lastrange != -1) {
      double logp = Math.log((double) (h-l) / lastrange) / Tools.LN2;
      System.err.println("\033[35mAC.stor "+l+", "+h
                        +" \t(+"+(-logp)+" bits)\033[m");
      if (h > lastrange) {
        warn("requested region exceeds last known range!  [diff="+(h-lastrange)+"]");
      }
    } else {
      System.err.println("\033[35mAC.stor "+l+", "+h+"\033[m");
    }
    coder.storeRegion(l,h);
    lastrange = -1; // discard last known range
  }

  public long getRange() {
    long res = coder.getRange();
    lastrange = res; // memorise last known range
    System.err.println("\033[33mAC.getR -> "+res+"\033[m");
    return res;
  }


  public void loadRegion(long l, long h, long t) {
    double logp = Math.log((double) (h-l) / t) / Tools.LN2;
    System.err.println("\033[32mAC.load "+l+", "+h+", "+t
                      +" \t\t("+logp+" bits)\033[m");
    coder.loadRegion(l,h,t);
    lastrange = -1; // discard last known range
  }

  public void loadRegion(long l, long h) {
    if (lastrange != -1) {
      double logp = Math.log((double) (h-l) / lastrange) / Tools.LN2;
      System.err.println("\033[32mAC.load "+l+", "+h
                        +" \t("+logp+" bits)\033[m");
      if (h > lastrange) {
        warn("requested region exceeds last known range!  [diff="+(h-lastrange)+"]");
      }
    } else {
      System.err.println("\033[32mAC.load "+l+", "+h+"\033[m");
    }
    coder.loadRegion(l,h);
    lastrange = -1; // discard last known range
  }

  public long getTarget(long t) {
    long res = coder.getTarget(t);
    System.err.println("\033[34mAC.getT "+t+" -> "+res+"\033[m");
    return res;
  }

  public long getTarget() {
    long res = coder.getTarget();
    System.err.println("\033[34mAC.getT -> "+res+"\033[m");
    return res;
  }


}
