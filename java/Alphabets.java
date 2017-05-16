
import java.util.Arrays;

/** A class with static methods for constructing various
  * uniform distributions over different alphabets. */
public class Alphabets {

  /** Constructs a uniform distribution over a selected
    * subset of printable ASCII letters and numbers,
    * so that each symbol takes exactly 5 bits to encode. */
  public static DiscreteUniform<Character> uniform5Bit() {
    Character[] syms = {
      '@',
      'A','B','C','D','E','F','G','H','I','J','K','L','M',
      'N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
      '[','\\',']','^','_'};
    return new DiscreteUniform<Character>(Arrays.asList(syms));
  }

  /** Constructs a uniform distribution over a selected
    * subset of printable ASCII letters and numbers,
    * so that each symbol takes exactly 6 bits to encode. */
  public static DiscreteUniform<Character> uniform6Bit() {
    Character[] syms = {
      ' ','0','1','2','3','4','5','6','7','8','9',
      'A','B','C','D','E','F','G','H','I','J','K','L','M',
      'N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
      'a','b','c','d','e','f','g','h','i','j','k','l','m',
      'n','o','p','q','r','s','t','u','v','w','x','y','z',
      '!' };
    return new DiscreteUniform<Character>(Arrays.asList(syms));
  }

  /** Constructs a uniform distribution over the 256
    * Unicode Braille characters. */
  public static DiscreteUniform<Character> uniformBraille() {
    Character[] syms = {
    '⠀','⠁','⠂','⠃','⠄','⠅','⠆','⠇','⠈','⠉','⠊','⠋','⠌','⠍','⠎','⠏',
    '⠐','⠑','⠒','⠓','⠔','⠕','⠖','⠗','⠘','⠙','⠚','⠛','⠜','⠝','⠞','⠟',
    '⠠','⠡','⠢','⠣','⠤','⠥','⠦','⠧','⠨','⠩','⠪','⠫','⠬','⠭','⠮','⠯',
    '⠰','⠱','⠲','⠳','⠴','⠵','⠶','⠷','⠸','⠹','⠺','⠻','⠼','⠽','⠾','⠿',
    '⡀','⡁','⡂','⡃','⡄','⡅','⡆','⡇','⡈','⡉','⡊','⡋','⡌','⡍','⡎','⡏',
    '⡐','⡑','⡒','⡓','⡔','⡕','⡖','⡗','⡘','⡙','⡚','⡛','⡜','⡝','⡞','⡟',
    '⡠','⡡','⡢','⡣','⡤','⡥','⡦','⡧','⡨','⡩','⡪','⡫','⡬','⡭','⡮','⡯',
    '⡰','⡱','⡲','⡳','⡴','⡵','⡶','⡷','⡸','⡹','⡺','⡻','⡼','⡽','⡾','⡿',
    '⢀','⢁','⢂','⢃','⢄','⢅','⢆','⢇','⢈','⢉','⢊','⢋','⢌','⢍','⢎','⢏',
    '⢐','⢑','⢒','⢓','⢔','⢕','⢖','⢗','⢘','⢙','⢚','⢛','⢜','⢝','⢞','⢟',
    '⢠','⢡','⢢','⢣','⢤','⢥','⢦','⢧','⢨','⢩','⢪','⢫','⢬','⢭','⢮','⢯',
    '⢰','⢱','⢲','⢳','⢴','⢵','⢶','⢷','⢸','⢹','⢺','⢻','⢼','⢽','⢾','⢿',
    '⣀','⣁','⣂','⣃','⣄','⣅','⣆','⣇','⣈','⣉','⣊','⣋','⣌','⣍','⣎','⣏',
    '⣐','⣑','⣒','⣓','⣔','⣕','⣖','⣗','⣘','⣙','⣚','⣛','⣜','⣝','⣞','⣟',
    '⣠','⣡','⣢','⣣','⣤','⣥','⣦','⣧','⣨','⣩','⣪','⣫','⣬','⣭','⣮','⣯',
    '⣰','⣱','⣲','⣳','⣴','⣵','⣶','⣷','⣸','⣹','⣺','⣻','⣼','⣽','⣾','⣿' };
    return new DiscreteUniform<Character>(Arrays.asList(syms));
  }

  /** Constructs a uniform distribution over 7-bit ASCII
    * characters.
    * @see UniformChar#ascii() */
  public static UniformChar uniform7Bit() {
    return UniformChar.ascii();
  }

  /** Constructs a uniform distribution over 256 byte values.
    * @see UniformByte */
  public static UniformByte uniform8Bit() {
    return new UniformByte();
  }

  /** Constructs a uniform distribution over the set of
    * ASCII digits (0123456789). */
  public static UniformChar uniformDigits() {
    return UniformChar.asciiDigits();
  }

}
