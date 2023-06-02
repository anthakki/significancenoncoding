
//
// A simple stream filter to transparently decompress
//  GZIP streams from files
//

public final class ZipFilter {

	private final static byte[] _GZIP_MAGIC = new byte[]{ (byte)0x1f, (byte)0x8b };
	private final static String _GZIP_SUFFIX = ".gz";

	// NB. compression level for .gz when writing
	public static int gzip_level = java.util.zip.Deflater.DEFAULT_COMPRESSION;

	private static boolean _peekMagic(java.io.BufferedInputStream stream, byte[] magic) throws java.io.IOException {
		stream.mark(magic.length);

		byte[] bytes = new byte[magic.length];
		boolean good = stream.read(bytes) == bytes.length;
		good = good && java.util.Arrays.equals(bytes, magic);

		stream.reset();

		return good;
	}

	public static java.io.InputStream filterInputStream(java.io.BufferedInputStream stream) throws java.io.IOException {
		if (_peekMagic(stream, _GZIP_MAGIC))
			return new java.util.zip.GZIPInputStream(stream);
		else
			return stream;
	}

	public static java.io.InputStream filterInputStream(java.io.InputStream stream) throws java.io.IOException {
		return filterInputStream(new java.io.BufferedInputStream(stream));
	}

	public static java.io.OutputStream filterOutputStream(java.io.OutputStream stream, java.io.File file) throws java.io.IOException {
		return filterOutputStream(stream, file.getPath());
	}

	public static java.io.OutputStream filterOutputStream(java.io.OutputStream stream, String name) throws java.io.IOException {
		if (name.endsWith(_GZIP_SUFFIX))
			return new java.util.zip.GZIPOutputStream(stream) {
				{
					this.def.setLevel(gzip_level);
				}
			};
		else
			return stream;
	}

} // ZipFilter
