/* $Id$ */

/** An interface for functions from elements of type A to elements
  * of type B. */
public interface Function<A,B> {

  public B eval(A a);

}
