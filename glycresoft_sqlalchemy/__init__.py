import warnings
from sqlalchemy import exc as sa_exc
warnings.simplefilter("ignore", category=sa_exc.SAWarning)
