# Quickstart

This tutorial shows how to fetch a paper with **PyPaperRetriever** and demonstrates writing clear docstrings.

## Fetch a paper

```python
from pypaperretriever import Retriever

retriever = Retriever()
paper = retriever.fetch("10.1038/nature12373")
print(paper.title)
```

## Docstring examples

```python
def calculate_sum(values):
    """Compute the sum of numeric values.

    Args:
        values (Iterable[float]): Numbers to be summed. An empty iterable returns ``0``.

    Returns:
        float: The total of all provided values.

    Raises:
        TypeError: If any element in ``values`` is not numeric.
    """
    return sum(values)


class DataProcessor:
    """Process and analyze numeric datasets.

    Attributes:
        scale (float): Factor applied to each element before aggregation.
    """

    def __init__(self, scale: float = 1.0):
        """Create a new processor.

        Args:
            scale (float, optional): Multiplier for data values. Defaults to ``1.0``.
        """
        self.scale = scale

    def process(self, data):
        """Scale and average a sequence of numbers.

        Args:
            data (Iterable[float]): Input values to process.

        Returns:
            float: Mean of the scaled values.

        Raises:
            ValueError: If ``data`` is empty.
        """
        if not data:
            raise ValueError("data must contain at least one element")
        scaled = [x * self.scale for x in data]
        return sum(scaled) / len(scaled)
```
