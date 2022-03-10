"""
Based on original code sourced from:
https://stackoverflow.com/questions/9727673/list-directory-tree-structure-in-python
- Original code Author: https://stackoverflow.com/users/2479038/abstrus

"""

from hashlib import md5
from pathlib import Path

class DisplayablePath(object):
    display_filename_prefix_middle = '├──'
    display_filename_prefix_last = '└──'
    display_parent_prefix_middle = '    '
    display_parent_prefix_last = '│   '

    def __init__(self, path, parent_path, is_last):
        self.path = Path(str(path))
        self.parent = parent_path
        self.is_last = is_last
        if self.parent:
            self.depth = self.parent.depth + 1
        else:
            self.depth = 0

    @property
    def displayname(self):
        if self.path.is_dir():
            return self.path.name + '/'
        return self.path.name

    @classmethod
    def make_tree(cls, root, parent=None, is_last=False, criteria=None):
        root = Path(str(root))
        criteria = criteria or cls._default_criteria

        displayable_root = cls(root, parent, is_last)
        yield displayable_root

        children = sorted(list(path
                               for path in root.iterdir()
                               if criteria(path)),
                          key=lambda s: str(s).lower())
        count = 1
        for path in children:
            is_last = count == len(children)
            if path.is_dir():
                yield from cls.make_tree(path,
                                         parent=displayable_root,
                                         is_last=is_last,
                                         criteria=criteria)
            else:
                yield cls(path, displayable_root, is_last)
            count += 1

    @classmethod
    def _default_criteria(cls, path):
        return True

    def _get_md5sum(self, path: Path):
        return md5(self.path.open('rb').read()).hexdigest()

    @property
    def displayname(self):
        if self.path.is_dir():
            return self.path.name + '/'
        return f"{self.path.name} md5: {self._get_md5sum(self.path)}"

    @property
    def md5sum(self):
        if self.path.is_dir():
            return None
        return self._get_md5sum(self.path)

    def displayable(self):
        if self.parent is None:
            return self.displayname

        _filename_prefix = (self.display_filename_prefix_last
                            if self.is_last
                            else self.display_filename_prefix_middle)

        parts = ['{!s} {!s}'.format(_filename_prefix,
                                    self.displayname)]

        parent = self.parent
        while parent and parent.parent is not None:
            parts.append(self.display_parent_prefix_middle
                         if parent.is_last
                         else self.display_parent_prefix_last)
            parent = parent.parent

        return ''.join(reversed(parts))

    @property
    def pandas_row(self):
        if self.path.is_dir():
            return None # don't return rows
        return {"filepath":str(self.path), "filename":str(self.path.name), "md5sum": self._get_md5sum(self.path)}
