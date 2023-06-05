
public abstract class ThreadEx extends java.lang.Thread {

	public void run() {
		try {
			runEx();
		}
		catch (java.lang.Exception error) {
			error.printStackTrace();
			System.exit(1);
		}
	}

	public abstract void runEx() throws java.lang.Exception;

}
