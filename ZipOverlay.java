
//
// A simple overlay to emulate reads from Zip files on
//  top of the filesystem for the reference data
//

public final class ZipOverlay {

	private static String _pathStringBasename(String path) {
		return new java.io.File(path).getName();
	}

	private static String _pathStringJoin(String path1, String path2) {
		return String.format("%s%s%s", path1, java.io.File.separator, path2);
	}

	private static String _pathStringStripTrailingSeparators(String path) {
		int len = path.length();
		for (; len > 0; --len)
			if (path.charAt(len - 1) != java.io.File.separatorChar)
				break;

		return path.substring(0, len);
	}

	private static interface _Object {
		public java.io.InputStream fileInputStream() throws java.io.IOException;

		public boolean exists();
	}

	private static class _FileObject implements _Object {
		protected java.io.File _file = null;

		public _FileObject(java.io.File file) {
			_file = file;
		}

		public java.io.InputStream fileInputStream() throws java.io.FileNotFoundException {
			return new java.io.FileInputStream(_file);
		}

		public boolean exists() {
			return _file.exists();
		}
	}

	private static class _ZipFileObject implements _Object {
		protected java.util.zip.ZipFile _zipFile = null;
		protected java.util.zip.ZipEntry _entry = null;

		public _ZipFileObject(java.util.zip.ZipFile zipFile, java.util.zip.ZipEntry entry) {
			_zipFile = zipFile;
			_entry = entry;
		}

		public java.io.InputStream fileInputStream() throws java.io.IOException {
			if (_entry.isDirectory())
				// NB. throw for directories like FileInputStream()
				throw new java.io.FileNotFoundException(String.format("%s (%s)", _pathStringBasename(_entry.getName()), "Is a directory"));

			return _zipFile.getInputStream(_entry);
		}

		public boolean exists() {
			return true;
		}
	}

	private static _Object _locateFile(java.io.File file) throws java.io.IOException {
		if (!file.exists()) {
			String tail = null;

			for (java.io.File parentFile = file;
					parentFile != null;
					parentFile = parentFile.getParentFile())
			{
				if (tail == null)
					tail = parentFile.getName();
				else
					tail = _pathStringJoin(parentFile.getName(), tail);

				java.io.File parentZipFile = new java.io.File(String.format("%s.zip", _pathStringStripTrailingSeparators(parentFile.getPath())));
				if (!parentZipFile.exists())
					continue;

				java.util.zip.ZipFile zipFile;
				try {
					zipFile = new java.util.zip.ZipFile(parentZipFile);
				}
				catch (java.util.zip.ZipException error) {
					// NB. not a zip or sth.. try next
					continue;
				}

				java.util.zip.ZipEntry entry = zipFile.getEntry(tail);
				if (entry == null)
					// NB. not in zip.. try next
					continue;

				return new _ZipFileObject(zipFile, entry);
			}
		}

		return new _FileObject(file);
	}

	public static java.io.InputStream fileInputStream(java.io.File file) throws java.io.IOException {
		return _locateFile(file).fileInputStream();
	}

	public static java.io.InputStream fileInputStream(String name) throws java.io.IOException {
		return fileInputStream(new java.io.File(name));
	}

	public static boolean exists(java.io.File file) {
		try {
			return _locateFile(file).exists();
		}
		catch (java.io.IOException error) {
			// NB. File.exists() doesn't throw either, not sure what it does..
			 // surely the result can be undetermined if i/o craps out
			throw new java.lang.RuntimeException(error);
		}
	}

	public static boolean exists(String name) {
		return exists(new java.io.File(name));
	}

} // ZipOverlay
