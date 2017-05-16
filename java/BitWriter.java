/* $Id$ */

import java.io.IOException;

/** An interface for writing low-level bit streams. */
public interface BitWriter {
  
  /** Writes a single bit (0 or 1).
    * The bit is encoded in the least-significant bit position
    * of the byte. */
  public void writeBit(byte bit) throws IOException;
  
  /** Writes a single bit.
    * @see Bit */
  public default void writeBit(Bit b) throws IOException {
    writeBit(b.byteValue());
  }
  
  /** Closes the bit stream.
    * Any buffered bits are flushed before closing. */
  public void close() throws IOException;

}
