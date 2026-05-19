# UppASD Portable Distribution

This is a portable distribution of UppASD that includes all necessary runtime libraries.

## Installation

### macOS
1. Open the DMG file
2. Drag UppASD to your Applications folder or any desired location
3. The application is ready to use

### From ZIP Archive
1. Extract the ZIP file to your desired location
2. The binary is in the `bin/` directory
3. Bundled libraries are in the `lib/` directory

## Running UppASD

### From Applications (macOS DMG)
- Double-click the UppASD application
- Or run from Terminal: `/Applications/UppASD.app/Contents/MacOS/uppasd`

### From Command Line
```bash
# Navigate to the extracted directory
cd /path/to/uppasd

# Run the binary
./bin/uppasd --help
```

## What's Included

This portable distribution includes:
- UppASD main binary (`uppasd`)
- Runtime libraries:
  - `libgfortran` (Fortran runtime)
  - `libquadmath` (Quadruple precision math)
  - `libgcc_s` (GCC runtime support)
  - `libgomp` (OpenMP runtime)

## System Requirements

- macOS 10.15 (Catalina) or later
- No additional software installation required
- All dependencies are bundled

## Verification

To verify the installation works correctly:

```bash
# Check dependencies (should show @rpath references)
otool -L ./bin/uppasd

# Test the binary
./bin/uppasd --version
```

## Troubleshooting

### "Cannot verify developer" error (macOS)
If you get a security warning:
1. Right-click the binary and select "Open"
2. Or go to System Preferences > Security & Privacy and click "Open Anyway"

### Library loading errors
If you see library loading errors:
1. Ensure the `lib/` directory is in the same parent directory as `bin/`
2. Check that library files have correct permissions: `chmod +x lib/*.dylib`

## Support

For issues and support, visit: https://github.com/UppASD/UppASD
