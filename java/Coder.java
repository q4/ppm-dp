/* $Id$ */

import java.io.IOException;

/** Interface for arithmetic coders. */
public interface Coder extends Encoder, Decoder {

    public void start_encode(BitWriter bw);
    public void finish_encode() throws IOException;

    public void start_decode(BitReader br) throws IOException;
    public void start_decode(BitReader br, boolean pad) throws IOException;
    public void finish_decode();

}
