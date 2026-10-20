"""Location of the input tables. Set RGBCR_DATA to override the default (<repository>/data)."""
import os

DATA = os.path.abspath(os.environ.get('RGBCR_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data'))) + os.sep
